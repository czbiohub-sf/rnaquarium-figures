#/bin/bash

# 02a_fixup_missing_fastp_price.sh
# fixup:
# older versions of the pipeline failed to output fastp and price stats
# which resulted in a csv frameshift
#
# cat data/intermediate/stats_merged_seqtype.csv | scripts/02a_fixup_missing_fastp_price.sh > data/intermediate/stats_merged_seqtype_fixup.csv


sed "s/^\([^,]*,[^,]*,[^,]*,[^,]*,[^,]*,\)\([^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*\),,,,,,$/\1,,,,,,\2/g" -