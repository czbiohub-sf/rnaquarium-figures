import os, math, pickle, argparse
# args
args = argparse.ArgumentParser()
args.add_argument("--input_dir", type=str, default="/path/to/Geneformer_RQfork/tokenize/RNAquarium/output")
args.add_argument("--loomfile", type=str, default="countsTable_60k_20240529.loom")
args.add_argument("--gtf_path", type=str, default="/path/to/Geneformer_RQfork/data/RNAquarium/GCF_049306965.1_GRCz12tu_genomic.gtf", help="full path to the gtf file")
args = args.parse_args()

data_dir = args.input_dir
input_file = os.path.join(data_dir, "tdigest_merged", "total_gene_median_dict.pickle")
gtf_path = args.gtf_path

outdir = os.path.join(data_dir, "token_dictionaries")
os.makedirs(outdir, exist_ok=True)

# read input 
with open(input_file, "rb") as fp:
    total_median_dict = pickle.load(fp)


########################
# ensembl_mapping_dict #
########################

ensembl_mapping_dict_RNAquarium = {k: k for k, v in total_median_dict.items() if not math.isnan(v)}
with open(os.path.join(outdir, "ensembl_mapping_dict_RNAquarium.pickle"), "wb") as fp:
    pickle.dump(ensembl_mapping_dict_RNAquarium, fp)


####################
# token_dictionary #
####################

token_dictionary = {'<pad>': 0,
                    '<mask>': 1}

for i, key in enumerate(total_median_dict, start=2):
    token_dictionary[key] = i

with open(os.path.join(outdir, "token_dictionary_RNAquarium.pickle"), "wb") as fp:
    pickle.dump(token_dictionary, fp)


#####################
# gene_name_id_dict #
#####################

def parse_gtf_attributes(attribute_string):
    """
    Parse the GTF attribute string and return a dictionary of key-value pairs.
    Example of attribute_string:
       'gene_id "ENSG00000198947"; transcript_id "ENST00000382353"; ...'
    """
    attr_dict = {}
    # Split on semicolon to get each key-value pair
    for attribute in attribute_string.strip().split(';'):
        attribute = attribute.strip()
        if not attribute:
            # Skip empty parts
            continue

        # Each attribute typically looks like key "value"
        # We can split on the first space to separate the key from the quoted value
        parts = attribute.split(' ', 1)
        if len(parts) != 2:
            continue

        key, value = parts
        # Remove surrounding quotes from the value
        value = value.strip('"')
        attr_dict[key] = value

    return attr_dict


def extract_gene_name(gtf_file):
    """
    Read a GTF file and return a dictionary mapping gene_id -> gene_biotype.
    """
    gene_name_dict = {}

    with open(gtf_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue

            # Split the GTF line into the 9 columns
            columns = line.strip().split('\t')
            if len(columns) < 9:
                continue

            feature_type = columns[2]  # e.g. "gene", "transcript", "exon", etc.
            attributes_str = columns[8]

            # We only want to look at lines describing a gene
            if feature_type != 'gene':
                continue

            # Parse the attributes into a dictionary
            attr_dict = parse_gtf_attributes(attributes_str)

            # Extract gene_id and gene_biotype if present
            if 'gene_id' in attr_dict:
                t_id = attr_dict['gene_id']
                # Some lines may not have gene_biotype, so use get() with a default
                gene_name = attr_dict.get('gene_name', 'N/A')
                gene_name_dict[gene_name] = t_id

    return gene_name_dict

gene_name_dict = extract_gene_name(gtf_path)

# save to pickle
with open(os.path.join(outdir, "gene_name_id_dict_RNAquarium.pickle"), "wb") as fp:
    pickle.dump(gene_name_dict, fp)