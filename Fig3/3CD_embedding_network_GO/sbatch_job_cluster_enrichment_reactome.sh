#!/bin/bash

#SBATCH --job-name=Geneformer_embeddings
#SBATCH --time=12:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH -e slurm.out/slurm-%A_%a.err
#SBATCH -o slurm.out/slurm-%A_%a.out

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

embeddings_dir="/path/to/Geneformer_RQfork/embeddings/extracted_embeddings"
layer=$1
python gp_enrich_cli.py \
  --embedding RQ="${embeddings_dir}/RNAquarium_${layer}.csv" \
  --embedding RQ_randomized="${embeddings_dir}/RNAquarium_randomized_${layer}.csv" \
  --embedding zebrahub-pretrained="${embeddings_dir}/RNAquarium_zebrahubData_${layer}.csv" \
  --embedding zebrahub-finetuned="${embeddings_dir}/RNAquarium_zebrahub_finetune_${layer}.csv" \
  --coexp-dir coexpression/leiden_results \
  --resolutions 0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10 \
  --cluster-spaces embedding \
  --go-sources REAC \
  --out-prefix "${layer}_Reactome" \
  --parallel \
  --max-workers 4 \
  --palette-map "RQ=#0a74b4,RQ_randomized=#707c92,zebrahub-pretrained=#9ba639,zebrahub-finetuned=#df6766"
# --make-ev-plot \
# --ev-max-pcs 255 \
# 