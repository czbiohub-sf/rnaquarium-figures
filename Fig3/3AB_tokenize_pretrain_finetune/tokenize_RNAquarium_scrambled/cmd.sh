module load mamba/23.1.0-3-pypy3
mamba activate Geneformer_RQfork

# define variables
script_dir=/path/to/Geneformer_RQfork/tokenize/scripts
counts_table_dir=/path/to/Geneformer_RQfork/data/RNAquarium
work_dir=/path/to/Geneformer_RQfork/tokenize/RNAquarium_scrambled
counts_table_basename=counts
gtf_basename=GCF_049306965.1_GRCz12tu_genomic
gene_id_to_ensembl_id_file="/path/to/Geneformer_RQfork/data/RNAquarium/xref/drerio_z12-geneid_to_z11-ENSDARG.csv"
output_prefix="RNAquarium_bulk"
additional_args="--scramble"
additional_prefix="_scrambled"

# derived variables
loom_file=$counts_table_basename.loom
gtf_path=$counts_table_dir/$gtf_basename.gtf
output_prefix="RNAquarium_bulk$additional_prefix"

# run the scripts
source "$script_dir/remove_output_and_logs.sh"
source "$script_dir/table_to_loom.sh"
source "$script_dir/shared_shell_commands.sh"
