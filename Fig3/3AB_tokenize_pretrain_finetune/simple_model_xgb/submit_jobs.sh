#!/bin/bash
base_path="/path/to/Geneformer_RQfork/simple_model/xgb/splits"
output_base_dir="/path/to/Geneformer_RQfork/simple_model/xgb/results"


train_eval_data_pair_count=0
train_data_not_exists=0
eval_data_not_exists=0
job_submitted=0
num_trials=50

dry_run=$1 # 1 for dry run, 0 for real run

# Iterate over all the pairs of train and eval data
# Outer loop: iterate from 0.05 to 0.95 with an increment of 0.05
for test_size in $(seq -f "%.2f" 0.05 0.05 0.95); do
    #echo "Current value: $test_size"
    #test_size=$(printf '%g' "$test_size")
    # Inner loop: iterate from 1 to 6
    for rep in {1..6}; do
        #echo "  Loop iteration: $i"
        for num_comp in 20 50 100; do
            #echo "test_size: $test_size, rep: $rep, num_comp: $num_comp"
             # Construct the filename based on the current loop values.
            eval_filename="eval_data_ncomps${num_comp}_testsize${test_size}_rep${rep}.pkl"
            eval_filepath="${base_path}/${eval_filename}"
      
            # Test if the file exists
            if [ ! -f "$eval_filepath" ]; then
                #echo "$eval_filepath does not exists."
                ((eval_data_not_exists++))
            fi

            train_filename="train_data_ncomps${num_comp}_testsize${test_size}_rep${rep}.pkl"
            train_filepath="${base_path}/${train_filename}"

            # Test if the file exists
            if [ ! -f "$train_filepath" ]; then
                echo "$train_filepath does not exists."
                ((train_data_not_exists++))
            fi
            ((train_eval_data_pair_count++))

            # Submit the job if both files exist
            if [ -f "$train_filepath" ] && [ -f "$eval_filepath" ] && [ $dry_run -eq 0 ]; then
                sbatch single_xgb_tune_and_train_job.sh $train_filepath $eval_filepath "${output_base_dir}/ncomps${num_comp}_testsize${test_size}_rep${rep}" $num_trials
                ((job_submitted++))
            fi
        done
    done
done


echo "Number of train-eval data pairs: $train_eval_data_pair_count"
echo "Number of train data that do not exist: $train_data_not_exists"
echo "Number of eval data that do not exist: $eval_data_not_exists"
echo "Number of jobs submitted: $job_submitted"