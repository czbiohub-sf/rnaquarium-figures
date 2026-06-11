import pickle, os, argparse, json
from utils import extract_gene_biotypes
import loompy
import pandas as pd
import numpy as np
import gzip

# args
args = argparse.ArgumentParser()
args.add_argument("--input_dir", type=str, default="/path/to/Geneformer_RQfork/data/RNAquarium", help="path to the input directory, contains the counts table and gtf file")
args.add_argument("--output_dir", type=str, default="./output")
args.add_argument("--counts_table_basename", type=str, default="countsTable_60k_20240529")
args.add_argument("--gtf_basename", type=str, default="GCF_049306965.1_GRCz12tu_genomic")
args.add_argument("--min_nonzero_threshold", type=int, default=7) # remove columns where the column has less than 7 non-zero values, Geneformer default is 7
args.add_argument("--scramble", action="store_true", help="scramble the counts table by shuffling the values") # scramble the counts table
args.add_argument("--randomize", action="store_true", help="randomize the counts table by replacing the values with random integer numbers") # randomize the counts table
args.add_argument("--gene_id_to_ensembl_id_file", type=str, default="/path/to/Geneformer_RQfork/data/RNAquarium/xref/drerio_z12-geneid_to_z11-ENSDARG.csv", help="path to the gene id to ensembl gene id mapping file")

# NEW: optional custom column attribute
args.add_argument("--attribute_name", type=str, default=None, help="name prefix for additional binary column attributes to add to col_attrs")
args.add_argument("--attribute_json", type=str, default=None, help="path to a JSON file mapping attribute labels to lists of column names")

args = args.parse_args()

# check args
if args.scramble and args.randomize:
    raise ValueError("Cannot scramble and randomize the counts table at the same time")

# directories
input_dir = args.input_dir
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

#####################
# load counts table #
#####################
counts_table_basename = args.counts_table_basename

# check if csv file exists
csv_file = os.path.join(input_dir, f"{counts_table_basename}.csv.gz")
if  os.path.exists(csv_file):
    # read the first line of a gzipped csv file and count the number of columns
    with gzip.open(csv_file, "rt") as f:
        first_line = f.readline()
        n_columns = len(first_line.split(","))
    # Specify data types for faster parsing
    dtype_mapping = {col: 'int32' for col in range(1, n_columns)}  # Adjust 'n_columns' to the number of columns in your data
    # load counts table
    print(f"reading csv file: {csv_file}", flush=True)
    counts_table = pd.read_csv(csv_file,
                                compression='gzip', 
                                index_col=0,            # Use the first column as row names
                                dtype=dtype_mapping,    # Ensure all other columns are integers
                                low_memory=False        # Disable internal type-checking chunks
                            )
    print(f"loaded the csv file, shape of counts table: {counts_table.shape}", flush=True)

# check if parquet file exists
parquet_file = os.path.join(input_dir, f"{counts_table_basename}.parquet")
if os.path.exists(parquet_file):
    # load counts table
    print(f"reading parquet file: {parquet_file}", flush=True)
    counts_table = pd.read_parquet(parquet_file)
    print(f"loaded the parquet file, shape of counts table: {counts_table.shape}", flush=True)

# set the first column as the index
counts_table.set_index(counts_table.columns[0], inplace=True)
counts_table.index.name = "Run"
counts_table.columns.name = "Gene"
print(counts_table.head(), flush=True)

if not os.path.exists(csv_file) and not os.path.exists(parquet_file):
    raise FileNotFoundError(f"Neither csv file nor parquet file found at path: {csv_file} or {parquet_file}")

# translate the gene ids to ensembl gene ids
gene_id_to_ensembl_id_df = pd.read_csv(args.gene_id_to_ensembl_id_file)
gene_id_to_ensembl_id_dict = dict(zip(gene_id_to_ensembl_id_df['gene_id'], gene_id_to_ensembl_id_df['Ensembl_gene_id']))
gene_id_to_ensembl_id_dict["Run"] = "Run" # preserve the run column
counts_table.columns = counts_table.columns.map(gene_id_to_ensembl_id_dict)
counts_table = counts_table.loc[:, counts_table.columns.notna()]
print(f"translated the gene ids to ensembl gene ids, shape of counts table: {counts_table.shape}", flush=True)

# transpose the counts table
print(counts_table.head(), flush=True)
counts_table = counts_table.T
print(f"transposed the counts table, shape of counts table: {counts_table.shape}", flush=True)
# print the first 5 rows and 5 columns of the counts table
print(counts_table.head(), flush=True)

##########################
# randomize counts table #
##########################
if args.randomize:
    print("Replace all values with random integers ...", flush=True)
    counts_table[:] = np.random.randint(0, 100, size=counts_table.shape, dtype=np.int32)

#########################
# scramble counts table #
#########################
if args.scramble:
    print("scrambling counts table...", flush=True)
    def shuffle_entire_dataframe(df):
        values = df.values.flatten() # flatten the dataframe into a 1D array
        n_shuffles = 3
        for _ in range(n_shuffles):
            np.random.shuffle(values) # shuffle the 1D array n_shuffles times
        return pd.DataFrame(values.reshape(df.shape), index=df.index, columns=df.columns) # reshape the 1D array into the original shape of the dataframe
        return pd.DataFrame(values.reshape(df.shape), index=df.index, columns=df.columns) # reshape the 1D array into the original shape of the dataframe

    counts_table = shuffle_entire_dataframe(counts_table)

########################
# filter counts table  #
########################
print("filtering counts table...", flush=True)
# remove columns where the column has less than 7 non-zero values
counts_table_filtered = counts_table.copy()
counts_table_filtered = counts_table_filtered.loc[:, (counts_table_filtered != 0).sum(axis=0) >= args.min_nonzero_threshold]
print(f"...removed columns (runs) where the column has less than {args.min_nonzero_threshold} non-zero values, shape of counts table: {counts_table_filtered.shape}", flush=True)

# drop rows where name starts with "_"
counts_table_filtered = counts_table_filtered[~counts_table_filtered.index.str.startswith('_')]
print(f"...dropped rows (genes/features) where name starts with '_', shape of counts table: {counts_table_filtered.shape}", flush=True)

##################
# load gene info #
##################
print("loading gene info...", flush=True)
gtf_path = os.path.join(input_dir, f"{args.gtf_basename}.gtf")
if not os.path.exists(gtf_path):
    raise FileNotFoundError(f"GTF file not found at path: {gtf_path}")
biotype_dict = extract_gene_biotypes(gtf_path)
print("done loading gene info", flush=True)
# translate the gene ids to ensembl gene ids
new_biotype_dict = {}
for key, value in biotype_dict.items():
    new_key = gene_id_to_ensembl_id_dict.get(key, key)
    new_biotype_dict[new_key] = value
biotype_dict = new_biotype_dict

########################################
# save the counts table in loom format #
########################################
print("saving the counts table in loom format...", flush=True)

# base col_attrs
col_attrs = {
    "cell": counts_table_filtered.columns.tolist(),
    "filter_pass": [1]*counts_table_filtered.shape[1], # set all cells to pass filter
    "n_counts": counts_table_filtered.sum(axis=0).tolist()
}

# NEW: add a binary attribute vector for each label in the JSON
# - JSON format: { "labelA": ["col1", "col3"], "labelB": ["col2"] }
# - For each label, we create a col_attrs key: f"{attribute_name}_{label}"
#   where value is a list of 1/0 aligned with counts_table_filtered.columns
if args.attribute_name is not None and args.attribute_json is not None:
    with open(args.attribute_json, "r") as f:
        label_to_columns = json.load(f)

    if not isinstance(label_to_columns, dict):
        raise ValueError("--attribute_json must point to a JSON object mapping labels to lists of column names")

    ordered_cols = counts_table_filtered.columns.tolist()

    for label, cols in label_to_columns.items():
        if not isinstance(cols, list):
            raise ValueError(f"Value for label '{label}' must be a list of column names")
        colset = set(cols)
        binary_vec = [1 if c in colset else 0 for c in ordered_cols]
        # name each attribute based on provided attribute_name prefix
        attr_key = f"{args.attribute_name}_{label.replace(' ', '_')}"
        col_attrs[attr_key] = binary_vec
        print(f"added attribute {attr_key} with length {len(binary_vec)}", flush=True)
        # add a string version of the attribute
        attr_key_str = f"{args.attribute_name}_{label.replace(' ', '_')}_str"
        col_attrs[attr_key_str] = [f"1" if c in colset else "0" for c in ordered_cols]
        print(f"added string attribute {attr_key_str} with length {len(col_attrs[attr_key_str])}", flush=True)

loompy.create(os.path.join(output_dir, f"{counts_table_basename}.loom"),
              np.asarray(counts_table_filtered),
              row_attrs={"ensembl_id": counts_table_filtered.index.tolist(),
                         "gene_type": [biotype_dict[gene] for gene in counts_table_filtered.index.tolist()]},
              col_attrs=col_attrs)
