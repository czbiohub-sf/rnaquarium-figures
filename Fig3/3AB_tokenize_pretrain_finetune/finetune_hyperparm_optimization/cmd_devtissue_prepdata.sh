#!/bin/bash

#SBATCH --job-name=Geneformer_dataprep
#SBATCH --time=1-00:00:00
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=8G
#SBATCH -e slurm.out/slurm-%A_%a.err
#SBATCH -o slurm.out/slurm-%A_%a.out

# load conda env
module load anaconda/latest
conda activate Geneformer_RQfork
module load cuda/12.4.0_550.54.14 # need to match cuda version of pytorch or will have GPU communication issues with nvcc
export CUDA_HOME=/hpc/apps/x86_64/cuda/12.4.0_550.54.14
module load cudnn/9.0.0.3_cuda12
module load hpcx/2.19

# source pretrained model paths
source /path/to/Geneformer_RQfork/pretrain/initial_weights/pretrained_models.sh
source /path/to/Geneformer_RQfork/pretrain/RNAquarium/pretrained_models.sh

# define model paths
initial_weights_mode_to_use="model_initial_weights_2048"
RQ_model_to_use="model_2048_E30"
RQ_model_to_use_randomized="model_randomized_2048_E30"
RQ_model_to_use_scrambled="model_scrambled_2048_E30"

# hardcode the model path for human weights
#model_human_weights="/path/to/Geneformer_RQfork/gf-6L-30M-i2048"

# use variable name expansion to get other model paths
model_initial_weights="/path/to/Geneformer_RQfork/pretrain/initial_weights/models/${!initial_weights_mode_to_use}/models"
model_RNAquarium="/path/to/Geneformer_RQfork/pretrain/RNAquarium/models/${!RQ_model_to_use}/models"
model_random="/path/to/Geneformer_RQfork/pretrain/RNAquarium/models/${!RQ_model_to_use_randomized}/models"
model_scrambled="/path/to/Geneformer_RQfork/pretrain/RNAquarium/models/${!RQ_model_to_use_scrambled}/models"

# define variables
# >>> check this path <<<
working_dir=/path/to/Geneformer_RQfork/finetune/input2048_30epochs_zebrahub_100-3000_tune # <<< check this
finetune_script="/path/to/Geneformer_RQfork/finetune/finetune.py"
freeze_layers=6
finetune_data_path="/path/to/Geneformer_RQfork/tokenize/zebrahub_100-3000/output/tokenized_data/zebrahub_devtissue_100-3000.dataset" # <<< check this
#no need to filter data anymore b/c the data is already filtered
#filter_data_dict='{"developmentalStage_anatomyOntologyClass":["larval-3dpf|mesenchyme","30 somites|central_nervous_system","larval-2dpf|central_nervous_system","larval-2dpf|mesenchyme","20 somites|central_nervous_system","05 somites|central_nervous_system","larval-5dpf|central_nervous_system","15 somites|central_nervous_system","larval-2dpf|periderm","larval-2dpf|lateral_mesoderm","30 somites|paraxial_mesoderm","larval-5dpf|paraxial_mesoderm","10 somites|central_nervous_system","larval-5dpf|intermediate_mesoderm","30 somites|lateral_mesoderm","larval-5dpf|hematopoietic_system","larval-10dpf|hematopoietic_system","05 somites|paraxial_mesoderm","larval-2dpf|paraxial_mesoderm","30 somites|periderm","larval-10dpf|periderm","larval-10dpf|paraxial_mesoderm","larval-3dpf|periderm","15 somites|periderm","30 somites|hematopoietic_system","20 somites|periderm","larval-10dpf|intermediate_mesoderm","15 somites|paraxial_mesoderm"]}'
#confu_plot_class_order="larval-3dpf|mesenchyme,30 somites|central_nervous_system,larval-2dpf|central_nervous_system,larval-2dpf|mesenchyme,20 somites|central_nervous_system,05 somites|central_nervous_system,larval-5dpf|central_nervous_system,15 somites|central_nervous_system,larval-2dpf|periderm,larval-2dpf|lateral_mesoderm,30 somites|paraxial_mesoderm,larval-5dpf|paraxial_mesoderm,10 somites|central_nervous_system,larval-5dpf|intermediate_mesoderm,30 somites|lateral_mesoderm,larval-5dpf|hematopoietic_system,larval-10dpf|hematopoietic_system,05 somites|paraxial_mesoderm,larval-2dpf|paraxial_mesoderm,30 somites|periderm,larval-10dpf|periderm,larval-10dpf|paraxial_mesoderm,larval-3dpf|periderm,15 somites|periderm,30 somites|hematopoietic_system,20 somites|periderm,larval-10dpf|intermediate_mesoderm,15 somites|paraxial_mesoderm"
cell_state_dict='{"state_key": "developmentalStage_anatomyOntologyClass", "states": "all"}'
finetune_task="devtissue"
attr_to_balance=""
attr_to_split="fish" # not used b/c will results in assertion error: assert len(split_attr_ids) == len(set(split_attr_ids))
n_cpu=12
fingerprint=$1

working_dir=$working_dir/scan_${finetune_task}
# # if exists, remove it
# if [ -d "$working_dir" ]; then
#     rm -rf $working_dir
# fi
mkdir -p $working_dir

# Loop through different test sizes
for test_size in $(seq 0.05 0.05 0.95); do
    perc=$(awk -v n="$test_size" 'BEGIN {printf "%.0fperc", n * 100}')
    
    echo "Running finetuning with test_size = $test_size ($perc)..."

    # finetune models
    cd $working_dir
    echo "#######################################################"
    echo "# Finetuning model pretrained on RNAquarium real data #"
    echo "#######################################################"
    python $finetune_script --working_dir . \
                        --output_prefix "zebrahub_${finetune_task}_classifier_freeze${freeze_layers}layers_testsize${perc}" \
                        --pretrained_model_path "${model_RNAquarium}" \
                        --finetune_data_path "${finetune_data_path}" \
                        --filter_data_dict "$filter_data_dict" \
                        --cell_state_dict "$cell_state_dict" \
                        --confu_plot_height 18 \
                        --confu_plot_width 18 \
                        --freeze_layers $freeze_layers \
                        --test_size $test_size \
                        --n_procs $n_cpu \
                        --prep_data_only \
                        --fingerprint "$fingerprint" 2>&1 | tee "$working_dir/zebrahub_${finetune_task}_classifier_freeze${freeze_layers}layers_testsize${perc}.log"

    # echo "##################################################"
    # echo "# Finetuning model pretrained on randomized data #"
    # echo "##################################################"
    # pretrain_variant="_randomized"
    # python $finetune_script --working_dir . \
    #                     --output_prefix "zebrahub_${finetune_task}_classifier${pretrain_variant}_freeze${freeze_layers}layers_testsize${perc}" \
    #                     --pretrained_model_path "$model_random" \
    #                     --finetune_data_path "${finetune_data_path}" \
    #                     --filter_data_dict "$filter_data_dict" \
    #                     --cell_state_dict "$cell_state_dict" \
    #                     --confu_plot_height 18 \
    #                     --confu_plot_width 18 \
    #                     --freeze_layers $freeze_layers \
    #                     --test_size $test_size \
    #                     --n_procs $n_cpu \
    #                     --prep_data_only \
    #                     --fingerprint "$fingerprint" 2>&1 | tee "$working_dir/zebrahub_${finetune_task}_classifier${pretrain_variant}_freeze${freeze_layers}layers_testsize${perc}.log"

    # echo "#################################################"
    # echo "# Finetuning model pretrained on scrambled data #"
    # echo "#################################################"
    # pretrain_variant="_scrambled"
    # python $finetune_script --working_dir . \
    #                     --output_prefix "zebrahub_${finetune_task}_classifier${pretrain_variant}_freeze${freeze_layers}layers_testsize${perc}" \
    #                     --pretrained_model_path "$model_scrambled" \
    #                     --finetune_data_path "${finetune_data_path}" \
    #                     --filter_data_dict "$filter_data_dict" \
    #                     --cell_state_dict "$cell_state_dict" \
    #                     --confu_plot_height 18 \
    #                     --confu_plot_width 18 \
    #                     --freeze_layers $freeze_layers \
    #                     --test_size $test_size \
    #                     --n_procs $n_cpu \
    #                     --prep_data_only \
    #                     --fingerprint "$fingerprint" 2>&1 | tee "$working_dir/zebrahub_${finetune_task}_classifier${pretrain_variant}_freeze${freeze_layers}layers_testsize${perc}.log"

    # echo "#################################################"
    # echo "# Finetuning model with human 30M-i2048 weights #"
    # echo "#################################################"
    # pretrain_variant="_human"
    # python $finetune_script --working_dir . \
    #                     --output_prefix "zebrahub_${finetune_task}_classifier${pretrain_variant}_freeze${freeze_layers}layers_testsize${perc}" \
    #                     --pretrained_model_path "$model_RNAquarium" \
    #                     --pretrained_model_path_weights_transfer_from "$model_human_weights" \
    #                     --finetune_data_path "${finetune_data_path}" \
    #                     --filter_data_dict "$filter_data_dict" \
    #                     --cell_state_dict "$cell_state_dict" \
    #                     --confu_plot_height 18 \
    #                     --confu_plot_width 18 \
    #                     --freeze_layers $freeze_layers \
    #                     --test_size $test_size \
    #                     --n_procs $n_cpu \
    #                     --prep_data_only \
    #                     --fingerprint "$fingerprint" 2>&1 | tee "$working_dir/zebrahub_${finetune_task}_classifier${pretrain_variant}_freeze${freeze_layers}layers_testsize${perc}.log"
                        
    # echo "###########################################################"
    # echo "# Finetuning model with initial weights (no pre-training) #"
    # echo "###########################################################"
    # pretrain_variant="_initialweights"
    # python $finetune_script --working_dir . \
    #                     --output_prefix "zebrahub_${finetune_task}_classifier${pretrain_variant}_freeze${freeze_layers}layers_testsize${perc}" \
    #                     --pretrained_model_path "$model_initial_weights" \
    #                     --finetune_data_path "${finetune_data_path}" \
    #                     --filter_data_dict "$filter_data_dict" \
    #                     --cell_state_dict "$cell_state_dict" \
    #                     --confu_plot_height 18 \
    #                     --confu_plot_width 18 \
    #                     --freeze_layers $freeze_layers \
    #                     --test_size $test_size \
    #                     --n_procs $n_cpu \
    #                     --prep_data_only \
    #                     --fingerprint "$fingerprint" 2>&1 | tee "$working_dir/zebrahub_${finetune_task}_classifier${pretrain_variant}_freeze${freeze_layers}layers_testsize${perc}.log"
done

echo "All iterations completed!"


 






















