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

################
# tokenization #
################

# tokenize RNAquarium data
cd $Geneformer_RQfork_dir/tokenize/RNAquarium
bash cmd.sh # can do in one of the login nodes

# tokenize RNAquarium randomized data
cd $Geneformer_RQfork_dir/tokenize/RNAquarium_randomized
bash cmd.sh

# tokenize RNAquarium scrambed data
cd $Geneformer_RQfork_dir/tokenize/RNAquarium_scrambled
bash cmd.sh # can do in one of the login nodes

# tokenize zebrahub data
cd $Geneformer_RQfork_dir/tokenize/zebrahub_100-3000
bash cmd.sh # can do in one of the login nodes

############
# pretrain #
############

# pretrain RNAquarium data
cd $Geneformer_RQfork_dir/pretrain/RNAquarium
sbatch pretrain_30epoch_RQ.sh
sbatch pretrain_30epoch_RQrandom.sh
sbatch pretrain_30epoch_RQscramble.sh
# NOTE:manual update pretrained_models.sh with the correct model prefix (which includes a unique timestamp)

# pretrain zebrahub data
cd $Geneformer_RQfork_dir/pretrain/zebrahub
bash pretrain_30epoch_ZH.sh
sbatch pretrain_30epoch_zebrahub.sh

###########################
# finetune - optimization #
###########################

# finetune RNAquarium data
cd $Geneformer_RQfork_dir/finetune/input2048_30epochs_zebrahub_100-3000_tune
# NOTE: manually submit (1) the prepdata job, (2) the finetune hyperopt job
# process hyperopt results
python $Geneformer_RQfork_dir/finetune/parse_hyperopt.py --input_dir $Geneformer_RQfork_dir/finetune/input2048_30epochs_zebrahub_100-3000_tune/slurm.out

# finetune zebrahub data
cd $Geneformer_RQfork_dir/finetune/zebrahub_model_zebrahub_100-3000_tune
# NOTE:manually submit (1) the prepdata job, (2) the finetune hyperopt job
# process hyperopt results
python $Geneformer_RQfork_dir/finetune/parse_hyperopt.py --input_dir $Geneformer_RQfork_dir/finetune/zebrahub_model_zebrahub_100-3000_tune/slurm.out


############
# finetune #
############

cd $Geneformer_RQfork_dir/finetune/input2048_30epochs_zebrahub_100-3000_hyperopted
# prep data
bash cmd_devtissue_prepdata.sh
# tune
bash cmd_devtissue_finetuneOnly.sh

################
# simple model #
################
cd /path/to/Geneformer_RQfork/simple_model/xgb
# generate the data split for xgboost
h5ad_path="/path/to/Geneformer_RQfork/data/zebrahub/zf_atlas_full_v4_release.developmentalStage_anatomyOntologyClass.100-3000.h5ad"
python prep_data.py --input_h5ad $h5ad_path --n_comps 20 --reps 6
python prep_data.py --input_h5ad $h5ad_path --n_comps 50 --reps 6
python prep_data.py --input_h5ad $h5ad_path --n_comps 100 --reps 6
# dry run
bash submit_jobs.sh 1 # check dry run results in stdout
# real run
bash submit_jobs.sh 0

#################
# plot results #
#################

# run jupyter notebook 
/path/to/Geneformer_RQfork/finetune/plot_scan_results.ipynb


