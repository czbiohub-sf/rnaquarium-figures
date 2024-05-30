# RNAquarium Pipeline: zebrafish filtering figures

This repository contains experimental scripts for analyzing and documenting the
RNAquarium pipeline as run on zebrafish RNA-seq data from NCBI SRA in April
2024. (As well as supporting details from mosquito and bat runs)

## Usage
```sh
git clone <url>
cd rnaquarium-pipeline-figures
```

To run the notebooks use the provided (conda)[https://docs.conda.io/en/latest/] environment:
```
conda env create -n rnaquarium-pipeline -f environment.yml
conda activate rnaquarium-pipeline
```

set up `data/` directory:
```
make data
```
follow the prompts.

For a distinct run, global settings are configured by first running the
`00_set_globals.ipynb` notebook.


## What's in this repo

[comment]: # {REPO_TREE}
```
.
├── environment-dev.yml
├── environment.yml
├── LICENSE
├── Makefile
├── notebooks
│   ├── 00_set_globals.ipynb
│   └── 01_massage_data.ipynb
├── pyproject.toml
├── README.md
└── scripts
    ├── 01_merge_stats.sh
    ├── 02_add_seqmethod_calls.py
    ├── 02a_fixup_missing_fastp_price.sh
    ├── 02b_fixup_old_reverse_sc_calls.sh
    ├── 03_make_stats_meta_table.py
    └── util
        └── update_readme_tree.sh
```

# License
This project is licensed under the BSD 3-Clause license - see the LICENSE file
for details.
