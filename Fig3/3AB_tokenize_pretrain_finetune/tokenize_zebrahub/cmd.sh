module load mamba/23.1.0-3-pypy3
mamba activate Geneformer_RQfork

# define variables
# data and script directories
script_dir=/path/to/Geneformer_RQfork/tokenize/scripts
RNAquarium_counts_dir=/path/to/Geneformer_RQfork/data/RNAquarium
RNAquarium_tokenize_dir=/path/to/Geneformer_RQfork/tokenize/RNAquarium
h5ad_dir=/path/to/Geneformer_RQfork/data/zebrahub
h5ad_basename=zf_atlas_full_v4_release.developmentalStage_anatomyOntologyClass.100-3000 # <<< check this
gtf_basename=GCF_049306965.1_GRCz12tu_genomic
gene_id_to_ensembl_id_file="/path/to/Geneformer_RQfork/data/RNAquarium/xref/drerio_z12-geneid_to_z11-ENSDARG.csv"

work_dir=/path/to/Geneformer_RQfork/tokenize/zebrahub_100-3000 # <<< check this
output_prefix="zebrahub_devtissue_100-3000" # <<< check this
additional_prefix=""

custom_attr_name_dict='''{"zebrafish_anatomy_ontology_class": "zebrafish_anatomy_ontology_class",
                            "developmental_stage": "developmental_stage",
                            "timepoint": "timepoint",
                            "fish": "fish",
                            "unique_cell_id": "unique_cell_id",
                            "developmentalStage_anatomyOntologyClass": "developmentalStage_anatomyOntologyClass"}'''

# *** fix the token dictionary path ***
token_dict_path="/path/to/Geneformer_RQfork/tokenize/RNAquarium/output/token_dictionaries/token_dictionary_RNAquarium.pickle"

# derived variables
loom_file=$h5ad_basename.loom
gtf_path=$RNAquarium_counts_dir/$gtf_basename.gtf
output_prefix="${output_prefix}${additional_prefix}"
# run the scripts
source "$script_dir/remove_output_and_logs.sh"
source "$script_dir/h5ad_to_loom.sh"
source "$script_dir/shared_shell_commands.sh"
