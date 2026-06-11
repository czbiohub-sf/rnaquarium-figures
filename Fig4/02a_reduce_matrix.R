# === PUBLIC RELEASE COPY ===
# Set <DATA_DIR> below to the folder holding the large Salmon matrix (see README).

#!/usr/bin/env Rscript
# =============================================================================
# 02a_reduce_matrix.R  --  reduce the full Salmon count matrix to curated contigs
# =============================================================================
# Stage 2a of the streamlined RNAquarium virus pipeline. This is the ONLY
# big-memory step and is meant to run as an HPC batch / interactive job, NOT in
# a notebook. It performs the matrix reduction in a large-memory session.
#
# Suggested launch:
#   srun --pty --cpus-per-task 12 --time 6:00:00 --mem 480G bash -l
#   module load r/4.4        # or your Global R module
#   Rscript 02a_reduce_matrix.R
# (Adjust --mem to your matrix size; the full matrix is ~74k rows x ~180k cols.)
#
# WHAT IT DOES
#   The full Salmon matrix columns are named with the SEPTEMBER contig names
#   (clusterLCA + cluster-size suffix from the Sept clustering run). The curated
#   contigs from 01 carry the OCTOBER names (slightly different cluster sizes /
#   CLUSTER labels). Here we reconcile this by joining on the run-invariant contig core id.
#   Here we instead join on the run-invariant contig CORE id (PRJ..._CONTIG_N,
#   the part before the first '|'), via sept_oct_crosswalk.tsv from 01. We then:
#     1. read the matrix header only (cheap),
#     2. find which Salmon columns map (by core id) to a curated October contig,
#     3. column-subset the big .tsv.gz to just those + the run_name column,
#     4. RENAME the kept columns from Sept names -> Oct names (so 02b sees the
#        same October naming as curated_contigs.tsv),
#     5. write the reduced matrix.
#
# INPUTS
#   - salmon_counts_matrix.tsv.gz        (full Salmon matrix; salmonpath)
#   - sept_oct_crosswalk.tsv          (from 01: contig_core, oct_full_name, is_target)
#   - curated_contigs.tsv             (from 01: the canonical curated October contigs)
# OUTPUT
#   - salmon_counts_fishassociated_reduced.tsv.gz   -> 02b   # -> 02b
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})
Sys.setenv("VROOM_CONNECTION_SIZE" = 50 * 1024 * 1024)  # 50 MB, for wide headers

# -- CONFIG -- edit paths here -----------------------------------------------
workingpath  <- getwd()
outpath0     <- "<DATA_DIR>"
salmonpath   <- "<DATA_DIR>/salmon_steps"

salmon_matrix_file <- "salmon_counts_matrix.tsv.gz"                 # in salmonpath
crosswalk_file     <- "sept_oct_crosswalk.tsv"                   # in workingpath (from 01)
curated_file       <- "curated_contigs.tsv"                      # in workingpath (from 01)
out_reduced_file   <- "salmon_counts_fishassociated_reduced.tsv.gz"  # -> 02b
# ----------------------------------------------------------------------------

core_id <- function(x) str_extract(x, "^[^|]+")   # part before first '|'

# --- 1. Load crosswalk + curated contigs (small) ----------------------------
crosswalk <- read_tsv(file.path(workingpath, crosswalk_file), show_col_types = FALSE)
curated   <- read_tsv(file.path(workingpath, curated_file),   show_col_types = FALSE)

# Oct names we want in the end (one per curated contig)
oct_names <- unique(curated$contig_withLCA_withcluster)
cat("Curated October contigs to recover:", length(oct_names), "\n")

# crosswalk gives contig_core <-> oct_full_name (+ is_target carve-out)
# Build core -> oct_name lookup (targets keep their full name as 'core')
cw <- crosswalk %>% distinct(contig_core, oct_full_name)

# --- 2. Read the Salmon matrix header only ----------------------------------
cat("Reading Salmon matrix header...\n")
hdr <- fread(file.path(salmonpath, salmon_matrix_file), nrows = 0, showProgress = TRUE)
all_cols <- names(hdr)
run_col  <- all_cols[1]                 # first column is the run/sample id
data_cols <- all_cols[-1]
cat("Total matrix columns:", length(all_cols), "(", length(data_cols), "contig columns )\n")

# --- 3. Match Salmon columns -> curated Oct contigs, by CORE id -------------
sept_tbl <- tibble(sept_name = data_cols,
                   contig_core = if_else(str_starts(data_cols, "NTtarget_"),
                                         data_cols, core_id(data_cols)))

matched <- sept_tbl %>%
  inner_join(cw, by = "contig_core") %>%       # attach the desired Oct name
  distinct(sept_name, .keep_all = TRUE)

cat("Salmon columns matched to a curated contig:", nrow(matched), "of",
    length(oct_names), "wanted\n")

# --- 3b. SPRIVIVIRUS SPECIAL CASE -------------------------------------------
# In the September Salmon run, the entire Sprivivirus cyprinus cluster (~1839
# contigs in the curated set) is represented by only the TWO NTtarget reference
# genomes -- the individual PRJ..._CONTIG_N Sprivivirus contigs have no column.
# (This is the "all sprivivirus replaced with 2 reference genomes" step from old
# the reference-genome step.) So those ~1839 curated contigs correctly have no Salmon column;
# their counts live entirely under the 2 NTtarget_ columns. We recover those 2
# columns here and label them with the curated Sprivivirus name so 02b can pool
# them. We assign a representative Oct name so the downstream cluster-collapse
# (which keys on the curated virus) picks them up.
sprivi_curated_name <- "Spring viraemia of carp virus"
sprivi_targets <- data_cols[str_detect(data_cols, "^NTtarget_") &
                            str_detect(data_cols, "Sprivivirus_cyprinus|Spring_viraemia_of_carp")]
cat("Sprivivirus NTtarget columns found in matrix:", length(sprivi_targets), "\n")

if (length(sprivi_targets) > 0) {
  sprivi_rows <- tibble(
    sept_name     = sprivi_targets,
    contig_core   = sprivi_targets,           # targets key on full name
    oct_full_name = sprivi_targets            # keep their own (Oct == Sept for targets)
  )
  matched <- bind_rows(matched, sprivi_rows) %>% distinct(sept_name, .keep_all = TRUE)
  cat("Added", length(sprivi_targets),
      "Sprivivirus target columns -> counts will collapse to:", sprivi_curated_name, "\n")
  cat("Total columns now matched:", nrow(matched), "\n")
}
# 02b note: when collapsing columns to curated viruses, the 2 NTtarget_ columns above
# must be mapped to clusterLCA_curated == "Spring viraemia of carp virus" (their column
# name's 3rd '|' field is CLUSTER1_Sprivivirus_cyprinus). The ~1839 individual Sprivivirus
# contigs have no column and correctly contribute 0 -- counts come solely from these 2.

# Report any wanted Oct contigs with NO Salmon column.
oct_recovered <- unique(matched$oct_full_name)
oct_missing   <- setdiff(oct_names, oct_recovered)
n_sprivi_missing <- sum(str_detect(oct_missing, "Sprivivirus_cyprinus"))
cat("Curated contigs with no direct Salmon column:", length(oct_missing),
    "( of which Sprivivirus:", n_sprivi_missing, ")\n")
cat("  EXPECTED: the ~1839 Sprivivirus contigs are counted via the 2 NTtarget columns (step 3b),\n")
cat("  so they correctly have no individual column. Any NON-Sprivivirus missing are worth a look:\n")
non_sprivi_missing <- oct_missing[!str_detect(oct_missing, "Sprivivirus_cyprinus")]
cat("  Non-Sprivivirus missing:", length(non_sprivi_missing), "\n")
if (length(non_sprivi_missing) > 0) {
  print(head(non_sprivi_missing, 20))
  writeLines(non_sprivi_missing, file.path(workingpath, "reduce_matrix_missing_contigs.txt"))
}

# --- 4. Column-subset the big matrix (read only the kept Sept columns) -------
keep_sept <- matched$sept_name
cat("Reading reduced matrix (", length(keep_sept), "contig columns )...\n")
df_subset <- read_tsv(file.path(salmonpath, salmon_matrix_file),
                      col_select = all_of(c(run_col, keep_sept)),
                      show_col_types = FALSE, progress = TRUE)
cat("Loaded:", nrow(df_subset), "rows x", ncol(df_subset), "cols\n")

# Guard: no readr auto-named columns
auto <- names(df_subset)[grepl("^\\.\\.\\.\\d+$", names(df_subset))]
if (length(auto) > 0) { cat("WARNING auto-named cols:\n"); print(auto) }

# --- 5. Rename kept columns Sept -> Oct names -------------------------------
# (so 02b joins cleanly against curated_contigs.tsv's October names)
rename_vec <- setNames(matched$oct_full_name, matched$sept_name)  # old=sept -> new=oct
present <- intersect(names(df_subset), names(rename_vec))
names(df_subset)[match(present, names(df_subset))] <- rename_vec[present]
cat("Renamed", length(present), "columns from September to October names\n")

# --- 6. Write reduced matrix ------------------------------------------------
out_path <- file.path(workingpath, out_reduced_file)
write_tsv(df_subset, out_path)
cat("Wrote", out_reduced_file, ":", nrow(df_subset), "x", ncol(df_subset), "\n")
cat("Done. -> feeds 02b_counts.ipynb\n")
