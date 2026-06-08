#!/bin/bash

#SBATCH --job-name=GeneFormer_pretrain_RQ_scramble
#SBATCH --time=2-10:00:00
#SBATCH --partition=gpu
#SBATCH --gpus=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=4
#SBATCH -e slurm.out/slurm-%A_%a.err
#SBATCH -o slurm.out/slurm-%A_%a.out

export MASTER_PORT=$((10000 + RANDOM % 10000))
# load conda env
module load anaconda/latest
conda activate Geneformer_RQfork
# module load cuda/12.5.1_555.42.06
# export CUDA_HOME=/hpc/apps/x86_64/cuda/12.5.1_555.42.06
module load cuda/12.4.0_550.54.14 # need to match cuda version of pytorch or will have GPU communication issues with nvcc
export CUDA_HOME=/hpc/apps/x86_64/cuda/12.4.0_550.54.14
module load cudnn/9.0.0.3_cuda12

#setting directories
working_dir="/path/to/Geneformer_RQfork/pretrain/RNAquarium"
suffix="_scrambled"
tokenize_dir=/path/to/Geneformer_RQfork/tokenize/RNAquarium$suffix


#main 
cd $working_dir
deepspeed --num_gpus=4 --num_nodes=1 \
         pretrain_geneformer_w_deepspeed.py --epochs 30 \
                                            --max_input_size 2048 \
                                            --output_prefix RNAquarium_bulk$suffix \
                                            --dataset_path "$tokenize_dir/output/tokenized_data/RNAquarium_bulk$suffix.dataset" \
                                            --lengths_path "$tokenize_dir/output/tokenized_data/RNAquarium_bulk${suffix}_lengths.pkl" \
                                            --token_dictionary_path "$tokenize_dir/output/token_dictionaries/token_dictionary_RNAquarium.pickle"


