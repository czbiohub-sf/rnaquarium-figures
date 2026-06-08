import scanpy as sc
import anndata as ad
import pandas as pd
import pickle, os, argparse
from utils import extract_gene_biotypes
import loompy
import numpy as np

# args
args = argparse.ArgumentParser()
args.add_argument("--input_dir", type=str, default="/path/to/Geneformer/data/zebrahub", help="path to the input directory, contains the h5ad file")
args.add_argument("--output_dir", type=str, default="./output")
args.add_argument("--h5ad_basename", type=str, default="zf_atlas_full_v1_release")
args.add_argument("--gtf_path", type=str, default="/path/to/Geneformer/data/RNAquarium/Danio_rerio.GRCz11.108.gtf", help="full path to the gtf file")
args.add_argument("--RNAquarium_tokenize_dir", type=str, default="/path/to/Geneformer/tokenize/RNAquarium", help="path to the RNAquarium tokenize directory, contains the ensembl_mapping_dict_RNAquarium.pickle file")
args.add_argument("--gene_id_to_ensembl_id_file", type=str, default="/path/to/Geneformer_RQfork/data/RNAquarium/xref/drerio_z12-geneid_to_z11-ENSDARG.csv", help="path to the gene id to ensembl gene id mapping file")

args = args.parse_args()



# directories
os.makedirs(args.output_dir, exist_ok=True)

# load h5ad file
print(f"loading h5ad file: {args.h5ad_basename}.h5ad")
h5ad_file = os.path.join(args.input_dir, f"{args.h5ad_basename}.h5ad")
adata = ad.read_h5ad(h5ad_file)

#############################################
# check gene ID consistency with RNAquarium #
#############################################

# print gene IDs to file
gene_ids = adata.var["gene_ids"]
with open(os.path.join(args.output_dir, "ZebraHub_gene_ids.txt"), "w") as f:
    for gene_id in gene_ids:
        f.write(f"{gene_id}\n")

# load RNAquarium gene IDs
RNAquarium_tokenize_dir = args.RNAquarium_tokenize_dir
ensembl_mapping_dict = pickle.load(open(os.path.join(RNAquarium_tokenize_dir,"output", "token_dictionaries", "token_dictionary_RNAquarium.pickle"), "rb"))

# load Danio_rerio.GRCz11.108.gtf
gene_biotypes = extract_gene_biotypes(args.gtf_path)
# translate the gene ids to ensembl gene ids
gene_id_to_ensembl_id_df = pd.read_csv(args.gene_id_to_ensembl_id_file)
gene_id_to_ensembl_id_dict = dict(zip(gene_id_to_ensembl_id_df['gene_id'], gene_id_to_ensembl_id_df['Ensembl_gene_id']))
new_biotype_dict = {}
for key, value in gene_biotypes.items():
    new_key = gene_id_to_ensembl_id_dict.get(key, key)
    new_biotype_dict[new_key] = value
gene_biotypes = new_biotype_dict

# check if all zebrahub gene IDs are in the RNAquarium gene IDs
# since RNAquarium gene IDs contain only protein-coding and miRNA genes, we skip gene IDs according to the
#  gene_biotypes dictionary
genes_not_in_RNAquarium = []
with open(os.path.join(args.output_dir, "ZebraHub_not_in_RNAquarium.txt"), "w") as f:
    for gene_id in gene_ids:
        if gene_id not in ensembl_mapping_dict and gene_id not in gene_biotypes:
            f.write(f"{gene_id}\n")
            genes_not_in_RNAquarium.append(gene_id)

# print the number of genes not in RNAquarium
print(f"Number of genes not in RNAquarium: {len(genes_not_in_RNAquarium)}")
if len(genes_not_in_RNAquarium) > 0:
    print(f"Genes not in RNAquarium: {genes_not_in_RNAquarium}")

# remove genes not in RNAquarium
print(f"Original number of genes: {adata.shape[1]}")
print(f"Original number of cells: {adata.shape[0]}")
print(f"Removing {len(genes_not_in_RNAquarium)} genes not in RNAquarium")
mask = ~adata.var["gene_ids"].isin(genes_not_in_RNAquarium)
adata = adata[:, mask]
print(f"After removing genes not in RNAquarium, number of genes: {adata.shape[1]}")
print(f"After removing genes not in RNAquarium, number of cells: {adata.shape[0]}")

###################
# convert to loom #
###################
#save the count table in loom format
print(f"Saving zebrahub counts in loom format")
# remove old loom file if it exists
if os.path.exists(os.path.join(args.output_dir, f"{args.h5ad_basename}.loom")):
    os.remove(os.path.join(args.output_dir, f"{args.h5ad_basename}.loom"))
loompy.create(os.path.join(args.output_dir, f"{args.h5ad_basename}.loom"),
              adata.X.T.toarray(), # transpose the count matrix so that cells are columns and genes are rows
              row_attrs={"ensembl_id": adata.var["gene_ids"].tolist(),
                         "gene_type": [gene_biotypes[gene] for gene in adata.var["gene_ids"].tolist()]},

              col_attrs={"cell": adata.obs.index.tolist(),
                         "unique_cell_id": adata.obs.index.tolist(),
                         "filter_pass": [1]*adata.shape[0], # set all cells to pass filter
                         "n_counts": adata.obs["total_counts"].tolist(),
                         "developmental_stage": adata.obs["developmental_stage"].tolist(),
                         "zebrafish_anatomy_ontology_class": adata.obs["zebrafish_anatomy_ontology_class"].tolist(),
                         "timepoint": adata.obs["timepoint"].tolist(),
                         "fish": adata.obs["fish"].tolist(),
                         "developmentalStage_anatomyOntologyClass": [a + "|" + b for a, b in zip(adata.obs["developmental_stage"].tolist(), adata.obs["zebrafish_anatomy_ontology_class"].tolist())] ,
                         })
