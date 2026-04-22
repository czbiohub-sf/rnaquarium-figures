# Figure 1 - Pipeline behind RNAquarium

 - Pipeline overview
 - Metrics for various key results from pipeline
   - number of datasets/projects processed
   - distributions for fraction of non-host reads
   - number of single cell vs bulk (could be a table instead?)
 - species generalizability (supplemental) 

## Inputs

| Metadata ID "rosetta" mapping | `/hpc/projects/balla_group/sra_experiments/SRA_metadata/dec2025_75k_submitteradded/zf_rosetta.tsv` |
| Full combined SRA metadata | `hpc/projects/balla_group/sra_experiments/SRA_metadata/dec2025_75k_submitteradded/all_zf_dates_devstage_tissue_tech_curated.tsv` |
| Read filtering stats for surviving runs | `/hpc/projects/balla_group/sra_experiments/versioned_zf_output/75k_unstable/host_mapping/stats-merged.csv` |
| Seq-Detective mate filter determinations | `/hpc/projects/balla_group/sra_experiments/versioned_zf_output/75k_unstable/host_mapping/reports/seq-detective-judgement-summary-all.txt` |
| Seq-Detective per-run metrics JSON | `/hpc/projects/balla_group/sra_experiments/versioned_zf_output/75k_unstable/host_mapping/download/<RUN_ID>/seq-detective-stats.json` |
| Nextflow process trace stats merged from two runs | `/hpc/projects/balla_group/sra_experiments/versioned_zf_output/75k_unstable/host_mapping/reports/trace-merged-dangerously.txt` |


copy input data files to local `data/` directory.

## Outputs
| notebook: input/output metric summary table  | `filtering.ipynb` |
| notebook: pie chart of sequence technologies by run count | `technologies.ipynb` |
| notebook: pie chart of sequencing datasets by metadata age | TODO |
| notebook: pie chart of sequencing datasets by tissue | TODO |
| notebook: contour plot of mapped reads mate 1 x mate 2 Seq-Detective metrics, contour per bioproject | TODO |
| notebook: contour plot of mapped reads x gene sparsity, separate contour by mate, per bioproject | TODO |

## TODO
 - Using `polars` for data processing
 - Using `great_tables` for tables
 - Using seaborn `kdeplot` for contour plots
 - Taking appropriate notes in `.claude/`

1. Find discrepancies between available data sources.
  - e.g. SRA runs with read filtering data but not Seq-Detective filtering determination, SRA runs
    with all steps in trace but no filtering stats.
  - Categorize and count runs by dropout point.
2. Write clean, notebook-like python scripts for generating the specified contour plots.  Use
   rosetta mapping to group run results into bioproject distributions.
