#!/usr/bin/env python3

# 02_add_seqmethod_calls.py
# the zebrafish pipeline run did *not* include sequence method (sc/bulk) calls
# in the output stats, but this information could be parsed out from work directories.
# this script combines the two sources into one table.
#
# python3 scripts/02_add_seqmethod_calls.py data/intermediate/stats_merged.csv data/rnaquarium_seq_type_calls_3.csv > data/intermediate/stats_merged_seqtype.csv # noqa: E501

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(
        description="merge SRA filtering stats csv with ncbi metadata tsv")
    parser.add_argument('stats', type=Path)
    parser.add_argument('method_calls', type=Path)
    args = parser.parse_args()

    stats_df = pd.read_csv(args.stats, header=0)
    types_df = pd.read_csv(args.method_calls, header=0,
                          names=["acc", "seq_type", "mate2_median"],
                          usecols=["acc", "seq_type"]
                          )

    joined_df = pd.merge(types_df, stats_df,
                         how='right', left_on='acc', right_on='id'
                        ).drop(columns='acc')
    print(joined_df.to_csv(index=False))


if __name__ == "__main__":
    main()
