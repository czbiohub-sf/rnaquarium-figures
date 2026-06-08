#!/bin/bash

#SBATCH --job-name=knn_leiden
#SBATCH --time=8:00:00
#SBATCH --array=1-25%25
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=96G
#SBATCH --cpus-per-task=12
#SBATCH -e slurm.out/slurm-%A_%a.err
#SBATCH -o slurm.out/slurm-%A_%a.out

resolution_array=(0.1 0.2 0.4 0.6 0.8 1 2 4 6 8 10 12 14 15 16 18 20 30 40 50 60 70 80 90 100)

# generate 0-based index based on array task ID
declare -x idx=$(( ${SLURM_ARRAY_TASK_ID} -1)) 

# Define a cleanup function to remove temporary files
cleanup() {
    echo "Cleaning up temporary files..."
    rm -rf "/tmp/${SLURM_JOB_ID}"
}
# Register the cleanup function to run on EXIT (which includes most failure cases)
trap cleanup EXIT SIGINT SIGTERM


# load conda env
module load anaconda/latest
conda activate Geneformer_RQfork

embeddings_dir="/path/to/Geneformer_RQfork/embeddings"
#setting directories
working_dir=$embeddings_dir
cd $working_dir

#main 
resolution=${resolution_array[$idx]}
echo "resolution: $resolution"

# python coexpression/leiden.py \
# --parquet coexpression/counts_tmm-norm.parquet \
# --axis samples_by_genes \
# --out coexpression/leiden_results/leiden_clusters_${resolution}.csv \
# --method pynndescent \
# --corr-type spearman \
# --min-corr 0.1 --n-neighbors 15 \
# --resolution $resolution \
# --n-jobs 12

python coexpression/leiden.py \
--parquet coexpression/counts_tmm-norm.parquet \
--axis samples_by_genes \
--out coexpression/leiden_results/leiden_clusters_${resolution}.csv \
--method pynndescent \
--corr-type robust --winsorize-q 0.01 \
--min-corr 0.1 --n-neighbors 15 \
--resolution $resolution \
--n-jobs 12
