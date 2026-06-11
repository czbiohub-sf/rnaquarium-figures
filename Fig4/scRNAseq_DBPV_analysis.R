# Example analysis using data from PRJNA667451

library(monocle3)
library(dplyr)
library(ggplot2)
# load data from cellranger output, then combine to one object
cds72 <- load_cellranger_data("~/Desktop/scRNAseq_DBPV/PRJNA667451/SRR12774372")
cds91 <- load_cellranger_data("~/Desktop/scRNAseq_DBPV/PRJNA667451/SRR12774391")
cds <- combine_cds(list(cds72, cds91))
remove(cds72, cds91)

# define mitochondria genes from GRCz12tu annotations
mt_genes <- c("trnF", "mt-rnr1", "trnV", "mt-rnr2", "trnL1", "ND1", "trnl",
              "trnQ", "trnM", "ND2", "trnW", "trnA", "trnN", "trnC", "trnY",
              "COX1", "trnS1", "trnD", "COX2", "trnK", "ATP8", "ATP6", "COX3",
              "trnG", "ND3", "trnR", "ND4L", "ND4", "trnH", "trnS2", "trnL2",
              "ND5", "ND6", "trnE", "CTYB", "trnT", "trnP")

# Calculate QC metrics
cds <- detect_genes(cds)

# Calculate mitochondrial percentage for each cell
# Match genes from mt list that are present in the data, check number
mt_genes_present <- intersect(mt_genes, rowData(cds)$id)
cat("Found", length(mt_genes_present), "mitochondrial genes in dataset\n")

mt_counts <- Matrix::colSums(exprs(cds)[rowData(cds)$id %in% mt_genes_present, ])
total_counts <- Matrix::colSums(exprs(cds))
colData(cds)$percent_mt <- (mt_counts / total_counts) * 100

# Calculate total UMI counts and number of genes detected per cell
colData(cds)$n_umi <- Matrix::colSums(exprs(cds))
colData(cds)$n_genes <- Matrix::colSums(exprs(cds) > 0)

# Visualize distributions
hist(colData(cds)$percent_mt, breaks = 50, main = "Mitochondrial %")
hist(log10(colData(cds)$n_umi), breaks = 50, main = "Log10 UMI counts")
hist(log10(colData(cds)$n_genes), breaks = 50, main = "Log10 genes detected")

# Filter cells
max_mito_percent <- 40
min_umi <- 500
min_genes <- 100
max_umi <- 50000

valid_cells <- colData(cds)$percent_mt < max_mito_percent &
  colData(cds)$n_umi > min_umi &
  colData(cds)$n_umi < max_umi &
  colData(cds)$n_genes > min_genes

cds_filtered <- cds[, valid_cells]

cat("Original cells:", ncol(cds), "\n")
cat("Filtered cells:", ncol(cds_filtered), "\n")
cat("Cells removed:", ncol(cds) - ncol(cds_filtered), "\n")

# Calculate PCA, look at fraction of variation explained by each
cds_filtered <- preprocess_cds(cds_filtered, num_dim = 75)
plot_pc_variance_explained(cds_filtered)

# Reduce dimensions for plotting UMAP, look for batch effects
cds_filtered <- reduce_dimension(cds_filtered)
plot_cells(cds_filtered)
plot_cells(cds_filtered, color_cells_by = "sample")

# Perform batch correction
cds_filteredb <- align_cds(cds_filtered, num_dim = 75, alignment_group = "sample")
cds_filteredb <- reduce_dimension(cds_filteredb)
plot_cells(cds_filteredb, color_cells_by = "sample")

# Cluster cells and inspect labeling
cds_filteredb <- cluster_cells(cds_filteredb, resolution = 1e-5)
plot_cells(cds_filteredb, group_label_size = 6, cell_size = 1)

plot_cells(cds_filteredb,
           genes = "nPicornavirus",
           label_cell_groups = FALSE,cell_size = 1.5, group_label_size = 5, scale_to_range = TRUE,
           show_trajectory_graph = FALSE)

# Find markers for each cluster
marker_test_res <- top_markers(cds_filteredb, group_cells_by = "cluster", genes_to_test_per_group = 75,verbose = T, speedglm.maxiter = 100)
write.csv(marker_test_res, file = "marker_genes.csv")

# Rename clusters based on cell type markers
colData(cds_filteredb)$assigned_cell_type <- as.character(clusters(cds_filteredb))
colData(cds_filteredb)$assigned_cell_type <- dplyr::recode(colData(cds_filteredb)$assigned_cell_type,
                                                           "1"="Macrophage",
                                                           "2"="RBC",
                                                           "3"="Fibroblast",
                                                           "4"="Endocardium (Ventricle)",
                                                           "5"="RBC",
                                                           "6"="Smooth muscle",
                                                           "7"="Apoptotic",
                                                           "8"="T cell",
                                                           "9"="Endothelial cells (lyve1)",
                                                           "10"="B cell",
                                                           "11"="Cardiomyocytes (Atrium)",
                                                           "12"="Neutrophil",
                                                           "13"="Endocardium (Atrium)",
                                                           "14"="Cardiomyocytes (dediff.)",
                                                           "15"="Endothelial cells (apnln)",
                                                           "16"="Valve Fibroblasts",
                                                           "17"="Monocyte",
                                                           "18"="Myelin cell",
                                                           "19"="Eosinophil",
                                                           "20"="Perivascular cell",
                                                           "21"="pDC",
                                                           "22"="ILC2")

plot_cells(cds_filteredb, group_cells_by="cluster", color_cells_by="assigned_cell_type", cell_size = 1, group_label_size = 5)

# Plot virus abundance
plot_cells(cds_filteredb,
           genes = "DBPV",
           label_cell_groups = TRUE,color_cells_by = "assigned_cell_type",cell_size = 1.5, group_label_size = 5, scale_to_range = TRUE,
           show_trajectory_graph = FALSE)+viridis::scale_color_viridis(na.value = "white")



