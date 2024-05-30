#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd

METADATA_COLUMNS = ["run_accession", "SRAStudy", "study_title", "Experiment", "Sample",
                    "BioSample", "Submission", "bioproject", "mixedsinglepairedflage",
                    "runsperbioproject", "bioproject_title", "bioproject_desc",
                    "experiment_title", "experiment_desc", "SRAStudy_abstract", "Sample_desc",
                    "earliestdate", "has_publication", "pmid", "pmid2", "pmid3", "url1",
                    "url2","url3", "url4", "colswithdata_count", "alejandro_devstage",
                    "alejandro_tissue", "scvsbulk_parsed1", "scvsbulk_parsed2",
                    "scvsbulk_parsed3", "scvsbulk_curated", "scvsbulk_notes", "spots",
                    "avgLength", "R1", "R2", "r1_r2", "lengthdiff", "lengthdiff2", "size_MB",
                    "experiment_alias", "LibraryName", "library_name", "LibraryStrategy",
                    "LibrarySelection", "LoadDate", "Organism", "Study_Pubmed_id", "sex",
                    "ebi_access_type", "organism part"
                    ]
KEEP_COLUMNS = ["run_accession", "bioproject", "earliestdate", "pmid",
                "alejandro_devstage", "scvsbulk_parsed3", "scvsbulk_curated",
                "scvsbulk_notes", "spots", "avgLength", "R1", "R2", "r1_r2",
                "lengthdiff", "lengthdiff2", "size_MB", "LibraryStrategy",
                "LibrarySelection", "LoadDate", "Organism", "sex", "organism part"]

def main() -> None:
    parser = argparse.ArgumentParser(
        description="merge SRA filtering stats csv with ncbi metadata tsv")
    parser.add_argument('stats', type=Path,
                        default=Path("stats_merged.csv"))
    parser.add_argument('metadata', type=Path,
                        default=Path("rnaquarium_metadata.tsv"))
    args = parser.parse_args()

    stats_df = pd.read_csv(args.stats, header=0)
    meta_df = pd.read_csv(args.metadata, header=0,
                          na_values=["Uncategorized"],
                          names=METADATA_COLUMNS,
                          usecols=KEEP_COLUMNS)

    return stats_df, meta_df


if __name__ == "__main__":
    main()
