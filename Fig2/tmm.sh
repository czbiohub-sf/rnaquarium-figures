# copy data
cp $PATH_TO_h5AD_FILE $PATH_TO_DATA_FOLDER/$h5AD_FILE_NAME

# tmm normalization of counts matrix, output is log2(TMM-CPM + prior_count)
counts_path= $PATH_TO_DATA_FOLDER/counts.parquet
output_path= $PATH_TO_DATA_FOLDER/counts_tmm-norm.parquet

rm -rf $output_path

python Fig2/scripts/tmm.py \
$counts_path \
$output_path \
--orientation auto --drop-last-n-cols 0 --output-scale logcpm \
--gene-id-to-ensembl-id-file $PATH_TO_DATA_FOLDER/drerio_z12-geneid_to_z11-ENSDARG.csv
