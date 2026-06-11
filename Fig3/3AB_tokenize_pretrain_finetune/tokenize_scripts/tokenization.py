from geneformer import TranscriptomeTokenizer
import os, pickle, argparse, json
from datasets import load_from_disk

# args
args = argparse.ArgumentParser()
args.add_argument("--input_dir", type=str, default="/path/to/Geneformer_RQfork/tokenize/RNAquarium/output")
args.add_argument("--loomfile", type=str, default="countsTable_60k_20240529.loom")
args.add_argument("--output_prefix", type=str, default="RNAquarium_bulk")
args.add_argument("--custom_attr_name_dict", type=str, default=None)
args.add_argument("--column_attributes_json_file", type=str, default=None, help="Path to the column attributes json file, this is used to get the column attributes from the json file")
args.add_argument("--token_dict_path", type=str, default=None)
args.add_argument("--file_format", type=str, default="loom")
args = args.parse_args()


# directories
rootdir = args.input_dir
outputdir = os.path.join(rootdir, "tokenized_data")
os.makedirs(outputdir, exist_ok=True)
data_directory = rootdir
dict_dir = os.path.join(rootdir, "token_dictionaries")

# file paths
loom_file = args.loomfile
GENE_MEDIAN_FILE = os.path.join(rootdir, "tdigest_merged", "detected_gene_median_dict.pickle")
TOKEN_DICTIONARY_FILE = os.path.join(dict_dir, "token_dictionary_RNAquarium.pickle")
ENSEMBL_MAPPING_FILE = os.path.join(dict_dir, "ensembl_mapping_dict_RNAquarium.pickle")
output_prefix = args.output_prefix

# custom attribute name dict
if args.custom_attr_name_dict is not None:
    custom_attr_name_dict = json.loads(args.custom_attr_name_dict)
else:
    custom_attr_name_dict = None

if args.column_attributes_json_file is not None:
    if custom_attr_name_dict is None:
        custom_attr_name_dict = {}

    with open(args.column_attributes_json_file, "r") as f:
        label_to_columns = json.load(f)

    if not isinstance(label_to_columns, dict):
        raise ValueError("--attribute_json must point to a JSON object mapping labels to lists of column names")

    for label, cols in label_to_columns.items():
        # compute the for each attribute (same way as in table_to_loom.py)
        attr_key = f"infection_{label.replace(' ', '_')}"
        print(f"added attribute {attr_key} during tokenization", flush=True)
        custom_attr_name_dict[attr_key] = attr_key

        # add a string version of the attribute
        attr_key_str = f"infection_{label.replace(' ', '_')}_str"
        print(f"added string attribute {attr_key_str} during tokenization", flush=True)
        custom_attr_name_dict[attr_key_str] = attr_key_str

# custom token dictionary
if args.token_dict_path is not None:
    TOKEN_DICTIONARY_FILE = args.token_dict_path # override the default token dictionary file if args.token_dict_path is provided
    print(f"Using custom token dictionary file: {TOKEN_DICTIONARY_FILE}")

#####################
# tokenize the data #
#####################

tk = TranscriptomeTokenizer(        
        custom_attr_name_dict=custom_attr_name_dict, # default is None
        nproc=1, # default is 1. If set to 4, an error will occur and it's related to multiprocessing.
        chunk_size=512, # default is 512
        model_input_size=2048, # 30M model needs to be 2048
        special_token=False, # 30M model needs to be False
        collapse_gene_ids=True, # default is True
        gene_median_file = GENE_MEDIAN_FILE, # default is to use that for the GeneFormer 90M model 
        token_dictionary_file = TOKEN_DICTIONARY_FILE, # default is to use that for the GeneFormer 90M model 
        gene_mapping_file = ENSEMBL_MAPPING_FILE # default is to use that for the GeneFormer 90M model 
)

tokenize = True
if tokenize == True:
        tk.tokenize_data(data_directory, 
                        outputdir, 
                        output_prefix, 
                        file_format=args.file_format)

print("Tokenization complete.")

######################################################
# extract the "length" field from the tokenized data #
######################################################

dataset_path = os.path.join(outputdir, f"{output_prefix}.dataset")

if os.path.exists(dataset_path):
    # 1. Load the dataset from the Arrow file
    print(f"Loading tokenized dataset from {dataset_path}")
    dataset = load_from_disk(dataset_path)

    # 2. Collect all lengths into a list
    print(f"Extracting lengths from the dataset")
    lengths = [dataset[i]["length"] for i in range(len(dataset))]

    # 3.Save the lengths to a pickle file
    with open(os.path.join(outputdir, f"{output_prefix}_lengths.pkl"), "wb") as f:
        pickle.dump(lengths, f)

    print(f"Extracted {len(lengths)} lengths and saved to {os.path.join(outputdir, f'{output_prefix}_lengths.pkl')}")
else:
    print(f"Dataset not found at {dataset_path}")
