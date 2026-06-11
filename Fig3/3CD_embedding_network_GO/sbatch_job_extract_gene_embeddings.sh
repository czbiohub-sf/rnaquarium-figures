#!/bin/bash

#SBATCH --job-name=Geneformer_embeddings
#SBATCH --time=1:00:00
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --constraint="h100_80|a100_80|h200"
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH -e slurm.out/slurm-%A_%a.err
#SBATCH -o slurm.out/slurm-%A_%a.out

# NOTE that an 80GB-memeory GPU is required for this job, as GPUs with 48GB memory will have insufficient memory errors

# Define a cleanup function to remove temporary files
cleanup() {
    echo "Cleaning up temporary files..."
    rm -rf "/tmp/${SLURM_JOB_ID}"
}
# Register the cleanup function to run on EXIT (which includes most failure cases)
trap cleanup EXIT SIGINT SIGTERM


# load conda env
module load anaconda/latest
conda activate Geneformer2025
module load cuda/12.4.0_550.54.14 # need to match cuda version of pytorch or will have GPU communication issues with nvcc
export CUDA_HOME=/hpc/apps/x86_64/cuda/12.4.0_550.54.14
module load cudnn/9.0.0.3_cuda12
module load hpcx/2.19


# model paths
RQ_model_to_use="model_2048_E30"
RQ_model_to_use_randomized="model_randomized_2048_E30"
RQ_model_to_use_scrambled="model_scrambled_2048_E30"
zebrahub_model_to_use="250810_230429_geneformer_Zebrahub_100-3000_L6_emb256_SL2048_E30_B12_LR0.001_LSlinear_WU10000_Oadamw_DS1"

# use variable name expansion to get other model paths
source /path/to/Geneformer_RQfork/pretrain/RNAquarium/pretrained_models.sh
model_RNAquarium="/path/to/Geneformer_RQfork/pretrain/RNAquarium/models/${!RQ_model_to_use}/models"
model_random="/path/to/Geneformer_RQfork/pretrain/RNAquarium/models/${!RQ_model_to_use_randomized}/models"
model_scrambled="/path/to/Geneformer_RQfork/pretrain/RNAquarium/models/${!RQ_model_to_use_scrambled}/models"
model_zebrahubData="/path/to/Geneformer_RQfork/pretrain/Zebrahub/models/${zebrahub_model_to_use}/models" # model pretrained on zebrahub data
model_zebrahub_finetune="/path/to/Geneformer_RQfork/finetune/input2048_30epochs_zebrahub_100-3000_hyperopted/scan_devtissue/zebrahub_devtissue_classifier_freeze6layers_testsize20perc/25Aug14t10pm1/25Aug14t10pm1_geneformer_cellClassifier_zebrahub_devtissue_classifier_freeze6layers_testsize20perc/ksplit1"

python extract_gene_embeddings.py --layer 0 --model_path_RQ $model_RNAquarium --model_path_random $model_random --model_path_scrambled $model_scrambled --model_path_zebrahub $model_zebrahubData --model_path_zebrahub_finetune $model_zebrahub_finetune # extract last layer embeddings
python extract_gene_embeddings.py --layer -1 --model_path_RQ $model_RNAquarium --model_path_random $model_random --model_path_scrambled $model_scrambled --model_path_zebrahub $model_zebrahubData --model_path_zebrahub_finetune $model_zebrahub_finetune # extract second-to-last layer embeddings