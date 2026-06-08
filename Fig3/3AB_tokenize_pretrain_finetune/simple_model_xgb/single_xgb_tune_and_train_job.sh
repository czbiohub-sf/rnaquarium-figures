#!/bin/bash

#SBATCH --job-name=XGBoost_tune&train
#SBATCH --time=12:00:00
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH -e slurm.out/slurm-%A_%a.err
#SBATCH -o slurm.out/slurm-%A_%a.out

module load anaconda/latest
module load cuda/12.5.1_555.42.06
module load cudnn/9.0.0.3_cuda12
export CUDA_HOME=/hpc/apps/x86_64/cuda/12.5.1_555.42.06
conda activate tune_xgb

mkdir -p "$3"

python tune_xgb.py --use_gpu --train_data $1 --eval_data $2 --output_dir $3 --num_trials $4 \
                   2>&1 | tee "${3}/xgb_tune_and_train_job.log"