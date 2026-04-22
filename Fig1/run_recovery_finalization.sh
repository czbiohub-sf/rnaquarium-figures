#!/bin/bash
#
# Master script to run all recovery finalization steps in sequence
# This will take a long time due to NFS file operations!
#

set -e  # Exit on error

echo "======================================================================="
echo "RECOVERY FINALIZATION - Starting at $(date)"
echo "======================================================================="
echo

# Step 1: Copy recovered download stats files
echo "STEP 1/4: Copying recovered download stats files..."
echo "This will take a long time due to NFS operations."
echo
python3 copy_recovered_download_stats.py
echo
echo "Step 1 complete at $(date)"
echo

# Step 2: Create augmented seq-detective TSV
echo "STEP 2/4: Creating augmented seq-detective judgement summary..."
echo
python3 create_augmented_seqdetective.py
echo
echo "Step 2 complete at $(date)"
echo

# Step 3: Create stats-with-dropouts-enhanced.csv
echo "STEP 3/4: Creating stats-with-dropouts-enhanced.csv..."
echo
python3 create_stats_with_dropouts_enhanced.py
echo
echo "Step 3 complete at $(date)"
echo

# Step 4: Create augmented host-filtering summary
echo "STEP 4/4: Creating augmented host-filtering summary..."
echo
python3 create_augmented_host_filtering_summary.py
echo
echo "Step 4 complete at $(date)"
echo

echo "======================================================================="
echo "ALL RECOVERY FINALIZATION STEPS COMPLETE at $(date)"
echo "======================================================================="
echo
echo "Output files created:"
echo "  - data/75k_unstable/seq-detective-judgement-summary-augmented.txt"
echo "  - data/75k_unstable/stats-with-dropouts-enhanced.csv"
echo "  - data/75k_unstable/host-filtering.summary.after-recovery.txt"
echo "  - data/75k_unstable/download_copy.log (copy operation log)"
echo
echo "These augmented files can now be used for generating pipeline figures."
