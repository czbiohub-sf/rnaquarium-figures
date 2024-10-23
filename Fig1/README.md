# Figure 1 - Pipeline behind RNAquarium

 - Pipeline overview
 - Metrics for various key results from pipeline
   - number of datasets/projects processed
   - distributions for fraction of non-host reads
   - number of single cell vs bulk (could be a table instead?)
 - species generalizability (supplemental) 

## Inputs
don't remember...

## Outputs


## Usage (probably outdated)

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
