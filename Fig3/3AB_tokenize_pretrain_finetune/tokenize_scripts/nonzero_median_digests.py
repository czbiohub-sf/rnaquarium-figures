import os, argparse
import numpy as np
import loompy as lp
import pandas as pd
import crick
import pickle
import math
from tqdm import tqdm

# args
args = argparse.ArgumentParser()
args.add_argument("--input_dir", type=str, default="/path/to/Geneformer/tokenize/RNAquarium/output")
args.add_argument("--loomfile", type=str, default="countsTable_60k_20240529.loom")
args = args.parse_args()

# check args
loompath = os.path.join(args.input_dir, args.loomfile)
if not os.path.exists(loompath):
    raise FileNotFoundError(f"Loom file not found at path: {loompath}")

rootdir = args.input_dir
output_file = args.loomfile.replace(".loom", ".gene_median_digest_dict.pickle")
outdir = os.path.join(rootdir, "tdigest")
outdir2 = os.path.join(rootdir, "tdigest_merged")

os.makedirs(outdir, exist_ok=True)
os.makedirs(outdir2, exist_ok=True)

databases = ["db"] # this is a dummy variable, it's useful in the human 30M dataset to merge tdigests from different databases
#####################
# process loom file #
#####################
with lp.connect(loompath) as data:
    # define coordinates of protein-coding or miRNA genes
    coding_miRNA_loc = np.where((data.ra.gene_type == "protein_coding") | (data.ra.gene_type == "miRNA"))[0]
    coding_miRNA_genes = data.ra["ensembl_id"][coding_miRNA_loc]
    
    # initiate tdigests
    median_digests = [crick.tdigest.TDigest() for _ in range(len(coding_miRNA_loc))]
    
    # initiate progress meters
    progress = tqdm(total=len(coding_miRNA_loc))
    last_view_row = 0
    progress.update(0)

    # print the number of cells and genes in the loom file
    print(f"Number of cells/columns: {data.shape[1]}")
    print(f"Number of genes/rows: {data.shape[0]}")
    
    for (ix, selection, view) in data.scan(items=coding_miRNA_loc, axis=0):
        # define coordinates of cells passing filter
        filter_passed_loc = np.where(view.ca.filter_pass == 1)[0]
        subview = view.view[:, filter_passed_loc]
        # normalize by total counts per cell and multiply by 10,000 to allocate bits to precision
        subview_norm_array = subview[:,:]/subview.ca.n_counts*10_000
         # if integer, convert to float to prevent error with filling with nan
        if np.issubdtype(subview_norm_array.dtype, np.integer):
            subview_norm_array = subview_norm_array.astype(np.float32)
        # mask zeroes from distribution tdigest by filling with nan
        nonzero_data = np.ma.masked_equal(subview_norm_array, 0.0).filled(np.nan)
        # update tdigests
        [median_digests[i+last_view_row].update(nonzero_data[i,:]) for i in range(nonzero_data.shape[0])]
        # update progress meters
        progress.update(view.shape[0])
        last_view_row = last_view_row + view.shape[0]
        
median_digest_dict = dict(zip(coding_miRNA_genes, median_digests))
with open(os.path.join(outdir, output_file), "wb") as fp:
    pickle.dump(median_digest_dict, fp)


################################################
# merge tdigests (only tdigest for RNAquarium) #
################################################

# merge new tdigests into total tdigest dict
# merge new tdigests into total tdigest dict
def merge_digest(dict_key_ensembl_id, dict_value_tdigest, new_tdigest_dict):
    new_gene_tdigest = new_tdigest_dict.get(dict_key_ensembl_id)
    if new_gene_tdigest is not None:
        dict_value_tdigest.merge(new_gene_tdigest)
        return dict_value_tdigest
    elif new_gene_tdigest is None:
        return dict_value_tdigest

# use tdigest1.merge(tdigest2) to merge tdigest1, tdigest2, ...tdigestn
# then, extract median by tdigest1.quantile(0.5)


# obtain gene list
#gene_info = pd.read_csv("/path/to/gene_info_table.csv", index_col=0)
#func_gene_list = [i for i in gene_info[(gene_info["gene_type"] == "protein_coding") | (gene_info["gene_type"] == "miRNA")]["ensembl_id"]]
func_gene_list = list(median_digest_dict.keys())

# initiate tdigests
median_digests = [crick.tdigest.TDigest() for _ in range(len(func_gene_list))]
total_tdigest_dict = dict(zip(func_gene_list, median_digests))

# merge tdigests
for current_database in databases:
    tdigest_dir = outdir
    for subdir, dirs, files in os.walk(tdigest_dir):	
        for file in files:
            if file.endswith(".gene_median_digest_dict.pickle"):
                print(f"merging {file}")
                with open(os.path.join(tdigest_dir, file), "rb") as fp:
                    tdigest_dict = pickle.load(fp)
                total_tdigest_dict = {k: merge_digest(k,v,tdigest_dict) for k, v in total_tdigest_dict.items()}

# save dict of merged tdigests
with open(os.path.join(outdir2, "total_gene_tdigest_dict.pickle"), "wb") as fp:
    pickle.dump(total_tdigest_dict, fp)

# extract medians and save dict
total_median_dict = {k: v.quantile(0.5) for k, v in total_tdigest_dict.items()}
with open(os.path.join(outdir2, "total_gene_median_dict.pickle"), "wb") as fp:
    pickle.dump(total_median_dict, fp)

# save dict of only detected genes' medians    
detected_median_dict = {k: v for k, v in total_median_dict.items() if not math.isnan(v)}
with open(os.path.join(outdir2, "detected_gene_median_dict.pickle"), "wb") as fp:
    pickle.dump(detected_median_dict, fp)