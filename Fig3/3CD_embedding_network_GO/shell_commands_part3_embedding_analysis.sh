##############################
#define paths and load conda #
##############################
Geneformer_RQfork_dir="/path/to/Geneformer_RQfork"
RQoutput_dir_path="/path/to/sra_experiments/versioned_zf_output/75k_unstable/host_mapping/host_counts"
GTF_name="GCF_049306965.1_GRCz12tu_genomic.gtf"
GTF_path="$RQoutput_dir_path/$GTF_name"

# load conda
module load mamba/23.1.0-3-pypy3
# create env
ENV_NAME="Geneformer_RQfork"
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "✅  Environment '$ENV_NAME' exists — activating it."
else
    echo "🆕  Environment '$ENV_NAME' not found — creating it."
    mamba create -n $ENV_NAME python=3.10.16
fi
conda activate $ENV_NAME

######################
# embedding analysis #
######################
cd /path/to/Geneformer_RQfork/embeddings

# extract embeddings real run 
# NOTE: check paths in both the sbatch_job_extract_gene_embeddings.sh file and the extract_gene_embeddings.py file
# NOTE: add --dry_run to the python calls in the sbatch_job_extract_gene_embeddings.sh file
sbatch sbatch_job_extract_gene_embeddings.sh

# prior to enrichement analysis, we perform co-expression analysis (this will be a baseline)
# tmm normalization [for co-expression analysis], output is log2(TMM-CPM + prior_count)
python coexpression/tmm.py \
/path/to/Geneformer_RQfork/data/RNAquarium/counts.parquet \
coexpression/counts_tmm-norm.parquet \
--orientation auto --drop-last-n-cols 0 --output-scale logcpm \
--gene-id-to-ensembl-id-file /path/to/Geneformer_RQfork/data/RNAquarium/xref/drerio_z12-geneid_to_z11-ENSDARG.csv
# NOTE: output is updated to be log2(TMM-CPM + prior_count) (2025-09-14)

# leiden clustering [for co-expression analysis]
mkdir -p coexpression/leiden_results
sbatch coexpression/batch_leiden.sh # started on 2025-09-15 12:15am

# write dense correlation matrix to file (for network analysis)
# Pearson correlation
python coexpression/leiden.py --parquet coexpression/counts_tmm-norm.parquet \
  --axis samples_by_genes \
  --method corr --corr-type pearson \
  --corr-out coexpression/corr_matrix_pearson/corr_matrix.npz

# Spearman (rank-based, outlier-resistant)
python coexpression/leiden.py --parquet coexpression/counts_tmm-norm.parquet \
  --axis samples_by_genes \
  --method corr --corr-type spearman \
  --corr-out coexpression/corr_matrix_spearman/corr_matrix.npz

# robust + light winsorization
python coexpression/leiden.py --parquet coexpression/counts_tmm-norm.parquet \
  --axis samples_by_genes \
  --method corr --corr-type robust --winsorize-q 0.001 \
  --corr-out coexpression/corr_matrix_robust_winsorize_0.001/corr_matrix.npz


# cluster embeddings and enrichment anlaysis [will also conduct enrichment analysis on co-expression cluster results and plot as a baseline]
#sbatch sbatch_job_cluster_enrichment.sh "last_layer"

sbatch sbatch_job_cluster_enrichment.sh "second_to_last_layer"
sbatch sbatch_job_cluster_enrichment_reactome.sh "second_to_last_layer"


# (replot)
python gp_enrich_cli.py \
--mock-summary-csv enrich_results/second_to_last_layer_GO_summary.csv \
--out-prefix "replot" \
--palette-map "RQ=#DF6766,RQ_randomized=#707C92,zebrahub-pretrained=#9ba639,zebrahub-finetuned=#df6766" \
--plot-series "RQ,RQ_randomized,co-expression" \
#--mock-umap-csv enrich_results/second_to_last_layer_GO_umap_embeddings.csv 

# produce input for cytoscape (KNN)
for k in 1 2 3 4; do
    python save_graph_for_cytoscape.py --top-k $k # NOTE: input paths are hardcoded in the script
done

# output full embedding cosine similarity matrix
emb_file="RNAquarium_second_to_last_layer.csv" 
python get_similarity_ranking.py --emb-sim-only --output-dir cosine_similarity_matrices --skip-coexpr \
--embeddings extracted_embeddings/$emb_file \
--emb-sim-output-file cosine_sim_mat_$emb_file

emb_file="RNAquarium_randomized_second_to_last_layer.csv" 
python get_similarity_ranking.py --emb-sim-only --output-dir cosine_similarity_matrices --skip-coexpr \
--embeddings extracted_embeddings/$emb_file \
--emb-sim-output-file cosine_sim_mat_$emb_file

emb_file="RNAquarium_zebrahub_finetune_second_to_last_layer.csv" 
python get_similarity_ranking.py --emb-sim-only --output-dir cosine_similarity_matrices --skip-coexpr \
--embeddings extracted_embeddings/$emb_file \
--emb-sim-output-file cosine_sim_mat_$emb_file

emb_file="esm2_15B_embeddings_gene_names.csv" 
esm_dir="/path/to/Geneformer_RQfork/embeddings/ESM/ESM_embeddings/esm2"
python get_similarity_ranking.py --emb-sim-only --output-dir cosine_similarity_matrices --skip-coexpr \
--embeddings $esm_dir/$emb_file \
--emb-sim-output-file cosine_sim_mat_$emb_file




# # get similarity rankings for selected genes
genelist=("ifnphi1" "stat1b" "usp18" "rsad2") 
mkdir -p gene_rankings
# embedding similarity
for gene in "${genelist[@]}"; do
    echo "Processing embedding similarity for: $gene"
    python get_similarity_ranking.py --gene $gene --skip-coexpr \
    --output-dir gene_rankings
done

# genelist=("isg15" "usp18" "ifnphi4" "ifng1" "ifng1r" "rsad2" "arg2" "stat1a" "stat1b" "stat2") 
# mkdir -p gene_rankings
# # embedding similarity
# for gene in "${genelist[@]}"; do
#     echo "Processing embedding similarity for: $gene"
#     python get_similarity_ranking.py --gene $gene --skip-coexpr \
#     --output-dir gene_rankings
# done

# genelist=("saa" "il1b" "il6" "tnfa") 
# mkdir -p gene_rankings
# # embedding similarity
# for gene in "${genelist[@]}"; do
#     echo "Processing embedding similarity for: $gene"
#     python get_similarity_ranking.py --gene $gene --skip-coexpr \
#     --output-dir gene_rankings
# done


# co-expression similarity (robust winsorization)
for gene in "${genelist[@]}"; do
  echo "Processing co-expression similarity (robust winsorization) for: $gene"
  python get_similarity_ranking.py --gene $gene --skip-embeddings \
    --coexpr coexpression/corr_matrix_robust_winsorize_0.001/corr_matrix.npz \
    --output-dir gene_rankings --output-prefix "robust_winsorize_0.001_"
  echo "Processing co-expression similarity (pearson) for: $gene"
  python get_similarity_ranking.py --gene $gene --skip-embeddings \
    --coexpr coexpression/corr_matrix_pearson/corr_matrix.npz \
    --output-dir gene_rankings --output-prefix "pearson_"
  echo "Processing co-expression similarity (spearman) for: $gene"
  python get_similarity_ranking.py --gene $gene --skip-embeddings \
    --coexpr coexpression/corr_matrix_spearman/corr_matrix.npz \
    --output-dir gene_rankings --output-prefix "spearman_"
done

# merge similarity rankings into a single table
for gene in "${genelist[@]}"; do
  echo "merging: $gene"
  python merge_similarity_ranking_table.py \
    "gene_rankings/${gene}_embedding_similarity_ranked.csv" \
    "gene_rankings/${gene}_robust_winsorize_0.001__coexpression_correlation_ranked.csv" \
    "gene_rankings/${gene}_pearson__coexpression_correlation_ranked.csv" \
    "gene_rankings/${gene}_spearman__coexpression_correlation_ranked.csv" \
    --output gene_rankings/merged_$gene.csv
done

# generate input for cytoscape
# take the top 50 neighbors for each gene and plot the connectivity network
emb_thresh_list=(0.43198 0.47143 0.47653, 0.44036)
coexpr_thresh_list=(0.57348 0.71827 0.6436 0.6655)

for i in "${!genelist[@]}"; do
  gene="${genelist[$i]}"
  emb_thresh="${emb_thresh_list[$i]}"
  coexpr_thresh="${coexpr_thresh_list[$i]}"
  echo "saving thresholded graph for: $gene with emb_thresh: $emb_thresh and coexpr_thresh: $coexpr_thresh"
  python save_thresholded_graph_for_cytoscape.py \
      --genes-csv "gene_rankings/merged_${gene}.csv" \
      --coexpr-npz  "/path/to/Geneformer_RQfork/embeddings/coexpression/corr_matrix_robust_winsorize_0.001/corr_matrix.npz"\
      --top-n 50 \
      --emb-thresh $emb_thresh \
      --coexpr-thresh $coexpr_thresh \
      --output-dir for_cytoscape \
      --output-prefix "graph_${gene}"
done

### produce multi-gene graph for cytoscape
# embedding threshold: 90 (retain only top 10% edges)
python save_multigene_graph_for_cytoscape.py \
    --genes ifnphi1 stat1b usp18 \
    --coexpr-npz  "/path/to/Geneformer_RQfork/embeddings/coexpression/corr_matrix_robust_winsorize_0.001/corr_matrix.npz"\
    --top-n 50 \
    --layers embeddings --emb-thresh 90 \
    --output-dir for_cytoscape \
    --output-prefix "graph_ifnphi1_stat1b_usp18"

python save_multigene_graph_for_cytoscape.py \
    --genes ifnphi1 rsad2 usp18 \
    --coexpr-npz  "/path/to/Geneformer_RQfork/embeddings/coexpression/corr_matrix_robust_winsorize_0.001/corr_matrix.npz"\
    --top-n 50 \
    --layers embeddings --emb-thresh 90 \
    --output-dir for_cytoscape \
    --output-prefix "graph_ifnphi1_rsad2_usp18"

python save_multigene_graph_for_cytoscape.py \
    --genes ifnphi1 rsad2 tnfa il10\
    --coexpr-npz  "/path/to/Geneformer_RQfork/embeddings/coexpression/corr_matrix_robust_winsorize_0.001/corr_matrix.npz"\
    --top-n 50 \
    --layers embeddings --emb-thresh 90 \
    --output-dir for_cytoscape \
    --output-prefix "graph_ifnphi1_rsad2_tnfa_il10"


python save_multigene_graph_for_cytoscape.py \
    --genes ifnphi1 rsad2 il6 il10 \
    --coexpr-npz  "/path/to/Geneformer_RQfork/embeddings/coexpression/corr_matrix_robust_winsorize_0.001/corr_matrix.npz"\
    --top-n 50 \
    --layers embeddings --emb-thresh 90 \
    --output-dir for_cytoscape \
    --output-prefix "graph_ifnphi1_rsad2_il6_il10"


# # get similarity rankings for mean vector of stat1 and stat2
# python get_similarity_ranking_composite.py --genes stat1a stat2 --skip-coexpr \
# --output-dir gene_rankings --output-prefix "stat1a_stat2_mean_" --embed-combine mean_l2


# # get similarity rankings for mean vector of stat1 and stat2
# python get_similarity_ranking_composite.py --genes stat1b stat2 --skip-coexpr \
# --output-dir gene_rankings --output-prefix "stat1b_stat2_mean_" --embed-combine mean_l2


##################
# NMP experiment #
##################
g0a=tbxta
g0b=sox2

python get_similarity_ranking_composite.py --genes "$g0a" "$g0b" --skip-coexpr \
--output-dir gene_rankings --output-prefix "${g0a}+${g0b}_mean_" --embed-combine mean_l2 # combine gene 1 and 2

# mesodermal track
g1=msgn1
g2=meox1
g3=myl1

python get_similarity_ranking_composite.py --genes "$g1" \
--query-vector-path "gene_rankings/${g0a}+${g0b}_mean__query_vector.csv" \
--query-vector-name "${g0a}+${g0b}" \
--embed-combine mean_l2 --skip-coexpr \
--output-dir gene_rankings --output-prefix "${g0a}+${g0b}+${g1}_mean_" # combine gene 0a 0b and 1

python get_similarity_ranking_composite.py --genes "$g2" \
--query-vector-path "gene_rankings/${g0a}+${g0b}+${g1}_mean__query_vector.csv" \
--query-vector-name "${g0a}+${g0b}+${g1}" \
--embed-combine mean_l2 --skip-coexpr \
--output-dir gene_rankings --output-prefix "${g0a}+${g0b}+${g1}+${g2}_mean_" # combine gene 0a 0b and 1 and 2

python get_similarity_ranking_composite.py --genes "$g3" \
--query-vector-path "gene_rankings/${g0a}+${g0b}+${g1}+${g2}_mean__query_vector.csv" \
--query-vector-name "${g0a}+${g0b}+${g1}+${g2}" \
--embed-combine mean_l2 --skip-coexpr \
--output-dir gene_rankings --output-prefix "${g0a}+${g0b}+${g1}+${g2}+${g3}_mean_" # combine gene 0a 0b and 1 and 2 and 3


# neural track
g1=pax6b
g2=pax6a

python get_similarity_ranking_composite.py --genes "$g1" \
--query-vector-path "gene_rankings/${g0a}+${g0b}_mean__query_vector.csv" \
--query-vector-name "${g0a}+${g0b}" \
--embed-combine mean_l2 --skip-coexpr \
--output-dir gene_rankings --output-prefix "${g0a}+${g0b}+${g1}_mean_" # combine gene 0a 0b and 1

python get_similarity_ranking_composite.py --genes "$g2" \
--query-vector-path "gene_rankings/${g0a}+${g0b}+${g1}_mean__query_vector.csv" \
--query-vector-name "${g0a}+${g0b}+${g1}" \
--embed-combine mean_l2 --skip-coexpr \
--output-dir gene_rankings --output-prefix "${g0a}+${g0b}+${g1}+${g2}_mean_" # combine gene 0a 0b and 1 and 2




####################################################
# NMP experiment - use tbxta as the starting point #
####################################################
# # mesodermal track
# g1=tbxta
# g2=msgn1
# g3=meox1
# g4=myl1

# python get_similarity_ranking_composite.py --genes "$g1"  --skip-coexpr \
# --output-dir gene_rankings 

# python get_similarity_ranking_composite.py --genes "$g2"  --skip-coexpr \
# --output-dir gene_rankings 

# python get_similarity_ranking_composite.py --genes "$g3"  --skip-coexpr \
# --output-dir gene_rankings 

# python get_similarity_ranking_composite.py --genes "$g4"  --skip-coexpr \
# --output-dir gene_rankings 

# python get_similarity_ranking_composite.py --genes "$g1" "$g2" --skip-coexpr \
# --output-dir gene_rankings --output-prefix "${g1}+${g2}_mean_" --embed-combine mean_l2 # combine gene 1 and 2

# python get_similarity_ranking_composite.py --genes "$g3" \
# --query-vector-path "gene_rankings/${g1}+${g2}_mean__query_vector.csv" \
# --query-vector-name "${g1}+${g2}" \
# --embed-combine mean_l2 --skip-coexpr \
# --output-dir gene_rankings --output-prefix "${g1}+${g2}+${g3}_mean_" # combine gene 1 and 2 and 3

# python get_similarity_ranking_composite.py --genes "$g4" \
# --query-vector-path "gene_rankings/${g1}+${g2}+${g3}_mean__query_vector.csv" \
# --query-vector-name "${g1}+${g2}+${g3}" \
# --embed-combine mean_l2 --skip-coexpr \
# --output-dir gene_rankings --output-prefix "${g1}+${g2}+${g3}+${g4}_mean_" # combine gene 1 and 2 and 3 and 4

# # neural track
# g1=tbxta
# g2=sox2
# g3=pax6b
# g4=pax6a

# python get_similarity_ranking_composite.py --genes "$g1"  --skip-coexpr \
# --output-dir gene_rankings 

# python get_similarity_ranking_composite.py --genes "$g2"  --skip-coexpr \
# --output-dir gene_rankings 

# python get_similarity_ranking_composite.py --genes "$g3"  --skip-coexpr \
# --output-dir gene_rankings 

# python get_similarity_ranking_composite.py --genes "$g4"  --skip-coexpr \
# --output-dir gene_rankings 

# python get_similarity_ranking_composite.py --genes "$g1" "$g2" --skip-coexpr \
# --output-dir gene_rankings --output-prefix "${g1}+${g2}_mean_" --embed-combine mean_l2 # combine gene 1 and 2

# python get_similarity_ranking_composite.py --genes "$g3" \
# --query-vector-path "gene_rankings/${g1}+${g2}_mean__query_vector.csv" \
# --query-vector-name "${g1}+${g2}" \
# --embed-combine mean_l2 --skip-coexpr \
# --output-dir gene_rankings --output-prefix "${g1}+${g2}+${g3}_mean_" # combine gene 1 and 2 and 3

# python get_similarity_ranking_composite.py --genes "$g4" \
# --query-vector-path "gene_rankings/${g1}+${g2}+${g3}_mean__query_vector.csv" \
# --query-vector-name "${g1}+${g2}+${g3}" \
# --embed-combine mean_l2 --skip-coexpr \
# --output-dir gene_rankings --output-prefix "${g1}+${g2}+${g3}+${g4}_mean_"

# # pairwise cosine similarity between mesodermal and neural tracks
# # loop over all composite vectors and compute against all composite vectors in the other track
# for g1 in tbxta tbxta+msgn1 tbxta+msgn1+meox1 tbxta+msgn1+meox1+myl1; do
#     echo "Processing: $g1"
#     python cosine.py "./gene_rankings/NMP_mesodermal_track/${g1}_mean__query_vector.csv" "./gene_rankings/NMP_neural_track/tbxta_mean__query_vector.csv"
#     python cosine.py "./gene_rankings/NMP_mesodermal_track/${g1}_mean__query_vector.csv" "./gene_rankings/NMP_neural_track/tbxta+sox2_mean__query_vector.csv"
#     python cosine.py "./gene_rankings/NMP_mesodermal_track/${g1}_mean__query_vector.csv" "./gene_rankings/NMP_neural_track/tbxta+sox2+pax6b_mean__query_vector.csv"
#     python cosine.py "./gene_rankings/NMP_mesodermal_track/${g1}_mean__query_vector.csv" "./gene_rankings/NMP_neural_track/tbxta+sox2+pax6b+pax6a_mean__query_vector.csv"
# done


# # get null distro for composite emb vectors
# python composite_nulls.py --mode random --combine sequential --n 4 --m 10000 \
#   --targets tbxta msgn1 meox1 myl1 sox2 pax6b pax6a --seed 123 \
#   --output-dir composite_nulls --output-prefix "sequential_" # sample 4 random genes and repeat 50 times, compute sequential composite emb vectors
# python composite_nulls.py --mode random --combine all --n 4 --m 10000 \
#   --targets tbxta msgn1 meox1 myl1 sox2 pax6b pax6a --seed 123 \
#   --output-dir composite_nulls --output-prefix "all_at_once_" # sample 4 random genes and repeat 50 times, compute all-at-once-mean composite emb vectors
# # plot null distro
# gene="msgn1"
# python plot_cosine_distribution.py composite_nulls/sequential__model_random_sequential_n4_m10000_cosine_to_targets.csv \
#   --columns msgn1 --bins 100 --density\
#   --mark 0.475:"tbxta+msgn1+meox1+myl1" --title "Null distribution of composite embedding cosine similarity to ${gene}"
# gene="meox1"
# python plot_cosine_distribution.py composite_nulls/sequential__model_random_sequential_n4_m10000_cosine_to_targets.csv \
#   --columns meox1 --bins 100 --density\
#   --mark 0.670:"tbxta+msgn1+meox1+myl1" --title "Null distribution of composite embedding cosine similarity to ${gene}"
# gene="myl1"
# python plot_cosine_distribution.py composite_nulls/sequential__model_random_sequential_n4_m10000_cosine_to_targets.csv \
#   --columns myl1 --bins 100 --density\
#   --mark 0.796:"tbxta+msgn1+meox1+myl1" --title "Null distribution of composite embedding cosine similarity to ${gene}"
# gene="sox2"
# python plot_cosine_distribution.py composite_nulls/sequential__model_random_sequential_n4_m10000_cosine_to_targets.csv \
#   --columns sox2 --bins 100 --density\
#   --mark 0.475:"tbxta+sox2+pax6b+pax6a" --title "Null distribution of composite embedding cosine similarity to ${gene}"
# gene="pax6b"
# python plot_cosine_distribution.py composite_nulls/sequential__model_random_sequential_n4_m10000_cosine_to_targets.csv \
#   --columns pax6b --bins 100 --density\
#   --mark 0.670:"tbxta+sox2+pax6b+pax6a" --title "Null distribution of composite embedding cosine similarity to ${gene}"
# gene="pax6a"
# python plot_cosine_distribution.py composite_nulls/sequential__model_random_sequential_n4_m10000_cosine_to_targets.csv \
#   --columns pax6a --bins 100 --density\
#   --mark 0.796:"tbxta+sox2+pax6b+pax6a" --title "Null distribution of composite embedding cosine similarity to ${gene}"

# # pairwise nulls
# python pairewise_nulls.py --n 2 --m 50000 --seed 123 --output-dir pairwise_nulls --output-prefix "pairwise_" 
# python plot_cosine_distribution.py pairwise_nulls/pairwise__pairwise_cosine_n2_m50000.csv \
#  --columns cosine_similarity --bins 100 --density --title "Null distribution of pairwise embedding cosine similarity"

