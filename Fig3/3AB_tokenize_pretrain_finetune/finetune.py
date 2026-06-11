# import torch.distributed as dist
# dist.init_process_group(backend='nccl') # this needs to be called before any other imports

import os
# os.environ["NCCL_DEBUG"] = "INFO" # this needs to be called before any other imports
# os.environ["OMPI_MCA_opal_cuda_support"] = "true"
# os.environ["CONDA_OVERRIDE_GLIBC"] = "2.56"

import datetime, json, time, argparse, gc
from geneformer import Classifier
import shutil
import pandas as pd
from datasets import load_from_disk

# args
args = argparse.ArgumentParser()
args.add_argument("--working_dir", type=str, default="./")
args.add_argument("--output_prefix", type=str, default="zebrahub_classifier_test")

args.add_argument("--pretrained_model_path", type=str, default="/path/to/Geneformer/pretrain/RNAquarium/models/250221_224411_geneformer_RNAquarium_bulk_L6_emb256_SL2048_E3_B12_LR0.001_LSlinear_WU10000_Oadamw_DS1/models", help="")
args.add_argument("--pretrained_model_path_weights_transfer_from", type=str, default=None, help="")
args.add_argument("--finetune_data_path", type=str, default="/path/to/Geneformer/tokenize/zebrahub/output/tokenized_data/zebrahub_sc.dataset")
args.add_argument("--eval_data_path", type=str, default=None, help="path to the evaluation data, need to set --special_eval, otherwise the test set (prepared by this script following --test_size) will be used for evaluation")

args.add_argument("--freeze_layers", type=int, help="number of layers to freeze", default=2)
args.add_argument("--test_size", type=float, help="size of the test set", default=0.2)

args.add_argument("--attr_to_balance", type=str, help="attribute to balance when splitting the data", default=None)
args.add_argument("--attr_to_split", type=str, help="attribute to split the data, e.g. 'patient_id' for splitting by patient while balancing other characteristics", default=None)
args.add_argument("--filter_data_dict", type=str, help="filters applied to the training data", default=None)
args.add_argument("--cell_state_dict", type=str, help="cell state dictionary", default=None)
args.add_argument("--max_samples_per_class", type=int, help="maximum number of samples per class", default=None)

args.add_argument("--confu_plot_class_order", type=str, default=None)
args.add_argument("--confu_plot_height", type=int, default=10)
args.add_argument("--confu_plot_width", type=int, default=10)

args.add_argument("--n_procs", type=int, default=8)
args.add_argument("--n_gpu", type=int, default=1)
args.add_argument("--use_tmp", action="store_true", help="use /tmp directory or use the directory set by --tmp_dir")
args.add_argument("--tmp_dir", type=str, default="/tmp", help="temporary directory")

args.add_argument("--local_rank", type=int, default=0, help="this is needed for deepspeed")

args.add_argument("--use_precomputed_hyperopt", action="store_true", help="use the precomputed hyperopt results")
args.add_argument("--hyperopt_results_path", type=str, default=None, help="path to the hyperopt results")
args.add_argument("--n_hyperopt_trials", type=int, default=0, help="number of hyperopt trials")
args.add_argument("--epochs", type=int, default=0.9, help="number of epochs for finetune")


args.add_argument("--prep_data_only", action="store_true", help="only prepare data, do not finetune model")
args.add_argument("--finetune_only", action="store_true", help="only finetune model, do not prepare data")
args.add_argument("--special_eval", action="store_true", help="evaluate model using the evaluation data provided by --eval_data_path")
args.add_argument("--fingerprint", type=str, default="", help="unique identifier for the run")
args.add_argument("--skip_evaluation", action="store_true", help="skip evaluation")
args.add_argument("--skip_confusion_plot", action="store_true", help="skip plotting confusion matrix")
args.add_argument("--force_finetune", action="store_true", help="force finetune even if model already exists")
args.add_argument("--force_eval", action="store_true", help="force evaluation even if evaluation results already exist")

args = args.parse_args()

# check arguments
if args.prep_data_only and args.finetune_only:
    raise ValueError("Cannot set both prep_data_only and finetune_only to True")

if any([args.prep_data_only, args.finetune_only]) and args.fingerprint == "":
    raise ValueError("fingerprint must be set if finetune_only or prep_data_only is True")

if args.special_eval and args.eval_data_path is None:
    raise ValueError("eval_data_path must be set if special_eval is True")

if args.use_precomputed_hyperopt and args.hyperopt_results_path is None:
    raise ValueError("hyperopt_results_path must be set if use_precomputed_hyperopt is True")

# set port
os.environ["MASTER_PORT"] = str(10000 + os.getpid() % 10000)

# directories
current_date = datetime.datetime.now()
datestamp = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}{current_date.hour:02d}{current_date.minute:02d}{current_date.second:02d}"
datestamp_min = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}"
if args.fingerprint != "":
    datestamp = args.fingerprint
    datestamp_min = args.fingerprint
working_dir = args.working_dir
output_dir = os.path.join(working_dir, f"{args.output_prefix}/{datestamp}")
os.makedirs(output_dir, exist_ok=True)
dataset_path = args.finetune_data_path

# functions
def get_directory_size(directory):
    """Returns the total size of a directory in bytes."""
    total_size = 0
    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            total_size += os.path.getsize(filepath)
    return total_size

def is_directory_above_threshold_mb(directory, threshold_mb):
    """Checks if the directory size exceeds a given threshold in MB."""
    threshold_bytes = threshold_mb * 1024**2  # Convert MB to Bytes
    dir_size = get_directory_size(directory)
    return dir_size > threshold_bytes, convert_size(dir_size)

def convert_size(size_bytes):
    """Converts bytes into a human-readable format (KB, MB, GB, etc.)."""
    if size_bytes == 0:
        return "0B"
    size_units = ["B", "KB", "MB", "GB", "TB"]
    i = 0
    while size_bytes >= 1024 and i < len(size_units) - 1:
        size_bytes /= 1024
        i += 1
    return f"{size_bytes:.2f} {size_units[i]}"

def read_csv_file(file_path):
    """
    Reads the CSV file at 'file_path' and returns a pandas DataFrame.
    """
    return pd.read_csv(file_path)

def find_row_by_value(df, target_column_name, target_column_value, extract_column_names):
    """
    Searches the DataFrame 'df' for rows where 'target_column_name' equals 'target_column_value'.
    Returns:
      - A dictionary of {column_name: row_value} for the first match
      - None if no match is found
    """
    df.rename(columns={"train_batch_size": "per_device_train_batch_size"}, inplace=True)
    matching_rows = df.loc[df[target_column_name] == target_column_value]
    if matching_rows.empty:
        return None  # or you can return an empty dict, or raise an exception
    # Extract the first matching row as a dictionary
    return matching_rows.iloc[0][extract_column_names].to_dict()

# read the hyperopt results
if args.use_precomputed_hyperopt:
    hyperopt_results_path = args.hyperopt_results_path # "/path/to/Geneformer/finetune/input2048_30epochs_zebrahub_200-2000_tune2_1-2-5epochs_200trials/hyperopt_results.csv"
    hyperopt_results = pd.read_csv(hyperopt_results_path)
    # find the test_size in the hyperopt results
    hyperparams = find_row_by_value(hyperopt_results, target_column_name = "test_size", target_column_value = args.test_size, 
                                    extract_column_names=["num_train_epochs", "learning_rate", "lr_scheduler_type", "warmup_steps", "weight_decay", "per_device_train_batch_size", "seed"])
    hyperparams["warmup_steps"] = int(hyperparams["warmup_steps"])
    hyperparams["seed"] = round(hyperparams["seed"])
    print(f"loaded hyperopt results for test_size {args.test_size}:\n {hyperparams}")
    # update the training args
    training_args = hyperparams
else:
    training_args = {
                    "num_train_epochs": 0.9,
                    "learning_rate": 0.000804,
                    "lr_scheduler_type": "polynomial",
                    "warmup_steps": 1812,
                    "weight_decay":0.258828,
                    "per_device_train_batch_size": 12,
                    "seed": 73,
                }

filter_data_dict=json.loads(args.filter_data_dict) if args.filter_data_dict else None
# OF NOTE: token_dictionary_file must be set to the gc-30M token dictionary if using a 30M series model
# (otherwise the Classifier will use the current default model dictionary)
# 30M token dictionary: https://huggingface.co/ctheodoris/Geneformer/blob/main/geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl
print("Initializing Geneformer Classifier", flush = True)
cc = Classifier(classifier="cell",
                cell_state_dict = json.loads(args.cell_state_dict) if args.cell_state_dict else None,
                filter_data=filter_data_dict,
                training_args=training_args,
                max_ncells=None,
                freeze_layers = args.freeze_layers,
                num_crossval_splits = 1,
                forward_batch_size=200,
                nproc=args.n_procs,
                ngpu=args.n_gpu)

###########################
# Prepare data            #
###########################
train_data_path = f"{output_dir}/{args.output_prefix}_labeled_train.dataset"
test_data_path = f"{output_dir}/{args.output_prefix}_labeled_test.dataset"
id_class_dict_path = f"{output_dir}/{args.output_prefix}_id_class_dict.pkl"
if not args.finetune_only: # skip preparing data if finetune only
    # check if prepared data exists
    exist_flag = all([os.path.exists(train_data_path), os.path.exists(test_data_path), os.path.exists(id_class_dict_path)])
    
    if exist_flag:
        size_flag = all([is_directory_above_threshold_mb(train_data_path, 10)[0], is_directory_above_threshold_mb(test_data_path, 10)[0], os.path.getsize(id_class_dict_path) > 100])
    else:
        size_flag = False
        
    if exist_flag and size_flag:
        print(f"Prepared data already exists:\n...{train_data_path}\n...{test_data_path}\n...{id_class_dict_path}", flush = True)
        print(f"...train_data_size: {convert_size(get_directory_size(train_data_path))}\n...test_data_size: {convert_size(get_directory_size(test_data_path))}\n...id_class_dict_size: {convert_size(os.path.getsize(id_class_dict_path))}", flush = True)
        print("...skipping data preparation", flush = True)
    else:
        print("Preparing data", flush = True)
        start_time = time.time()
        cc.prepare_data(input_data_file=dataset_path,
                        output_directory=output_dir,
                        output_prefix=args.output_prefix,
                        test_size = args.test_size, # hold out data for testing
                        attr_to_balance=args.attr_to_balance, #List of attribute keys on which to balance data while splitting on attr_to_split
                        #attr_to_split=args.attr_to_split # *need to implement fish*,  Key for attribute on which to split data while balancing potential confounders e.g. "patient_id" for splitting by patient while balancing other characteristics
                        max_samples_per_class=args.max_samples_per_class
                    )
        end_time = time.time()
        print(f"Time taken to prepare data: {end_time - start_time} seconds", flush = True)
else:
    print("Skipping data preparation: finetune_only flag was set", flush = True)
    # check if prepared data exists
    print(f"checking prepared data", flush = True)
    if not os.path.exists(train_data_path):
        raise FileNotFoundError(f"File {train_data_path} does not exist")
    else:
        print(f"...prepared data file exists: {train_data_path}", flush = True)
    if not os.path.exists(id_class_dict_path):
        raise FileNotFoundError(f"File {id_class_dict_path} does not exist")
    else:
        print(f"...prepared id_class_dict file exists: {id_class_dict_path}", flush = True)
    if not os.path.exists(test_data_path):
        raise FileNotFoundError(f"File {test_data_path} does not exist")
    else:
        print(f"...prepared labeled test file exists: {test_data_path}", flush = True)

###########################
# Load model and finetune #
###########################
if not args.prep_data_only: # skip finetuning if prep_data_only
    
    model_dir_path = f"{output_dir}/{datestamp_min}_geneformer_cellClassifier_{args.output_prefix}/ksplit1/"
    model_file_path = os.path.join(model_dir_path, "model.safetensors")

    # delete the finetuned model directory if force_finetune is set
    if args.force_finetune:
        print(f"force_finetune flag was set, deleting model directory: {os.path.dirname(os.path.dirname(model_dir_path))}", flush = True)
        if os.path.exists(os.path.dirname(os.path.dirname(model_dir_path))):
            shutil.rmtree(os.path.dirname(os.path.dirname(model_dir_path)))

    # check if finetuned model already exists
    if os.path.exists(model_file_path) and os.path.getsize(model_file_path) > 0:
        print(f"Finetuned model already exists: {model_file_path}", flush = True)
        print(f"...Finetuned model file size: {os.path.getsize(model_file_path) / (1024 * 1024):.2f} MB", flush=True)
        print("...skipping finetune training and jumping to evaluation", flush = True)
    else:
        print("finetuning model", flush = True)
        # remove the model directory if it exists
        if os.path.exists(os.path.dirname(model_dir_path)):
            shutil.rmtree(os.path.dirname(model_dir_path))
            shutil.rmtree(os.path.dirname(os.path.dirname(model_dir_path)))
        custom_stamp = args.fingerprint if args.fingerprint != "" else None

        if args.use_tmp: # if use_tmp is set, use the /tmp directory to store a copy of the input data for faster I/O
            # copy the input data to the /tmp directory
            print(f"using {args.tmp_dir} directory to store input/training data", flush = True)
            tmp_dir = os.path.join(args.tmp_dir, datestamp) # store the input data in a temporary directory with the datestamp
            if os.path.exists(os.path.join(tmp_dir, f"{args.output_prefix}_labeled_train.dataset")):
                shutil.rmtree(os.path.join(tmp_dir, f"{args.output_prefix}_labeled_train.dataset"))
            print(f"copying {os.path.join(output_dir, f'{args.output_prefix}_labeled_train.dataset')} to {os.path.join(tmp_dir, f'{args.output_prefix}_labeled_train.dataset')}", flush = True)
            shutil.copytree(f"{output_dir}/{args.output_prefix}_labeled_train.dataset", os.path.join(tmp_dir, f"{args.output_prefix}_labeled_train.dataset"))
            if os.path.exists(os.path.join(tmp_dir, f"{args.output_prefix}_id_class_dict.pkl")):
                os.remove(os.path.join(tmp_dir, f"{args.output_prefix}_id_class_dict.pkl"))
            print(f"copying {os.path.join(output_dir, f'{args.output_prefix}_id_class_dict.pkl')} to {os.path.join(tmp_dir, f'{args.output_prefix}_id_class_dict.pkl')}", flush = True)
            shutil.copy(f"{output_dir}/{args.output_prefix}_id_class_dict.pkl", os.path.join(tmp_dir, f"{args.output_prefix}_id_class_dict.pkl"))
            output_dir = tmp_dir
            # remove finetuned model directory if it exists
            if os.path.exists(os.path.join(tmp_dir, f"{datestamp_min}_geneformer_cellClassifier_{args.output_prefix}")):
                shutil.rmtree(os.path.join(tmp_dir, f"{datestamp_min}_geneformer_cellClassifier_{args.output_prefix}"))

        # print out the dimensions of the training data and the number of classes
        _data = load_from_disk(f"{output_dir}/{args.output_prefix}_labeled_train.dataset")
        print(f"Training data rows: {_data.num_rows}, columns: {len(_data.column_names)}", flush=True)
        print(f"Number of classes: {len(_data.unique('label'))}", flush=True)
        del _data
        gc.collect()

        start_time = time.time()
        all_metrics = cc.validate(model_directory=args.pretrained_model_path,
                                prepared_input_data_file=f"{output_dir}/{args.output_prefix}_labeled_train.dataset",
                                id_class_dict_file=f"{output_dir}/{args.output_prefix}_id_class_dict.pkl",
                                output_directory=output_dir,
                                output_prefix=args.output_prefix,
                                model_directory_weights_transfer_from=args.pretrained_model_path_weights_transfer_from,
                                custom_stamp=custom_stamp,
                                n_hyperopt_trials=args.n_hyperopt_trials
                                )
                                # to optimize hyperparameters, set n_hyperopt_trials=100 (or alternative desired # of trials)
        end_time = time.time()
        print(f"Time taken to finetune model: {end_time - start_time} seconds", flush = True)

        # copy the finetuned model to the output directory
        if args.use_tmp:
            output_dir = os.path.join(working_dir, f"{args.output_prefix}/{datestamp}/") # reset the output directory to the original output directory
            # if target directory exists, remove it
            if os.path.exists(os.path.join(output_dir, f"{datestamp_min}_geneformer_cellClassifier_{args.output_prefix}")):
                shutil.rmtree(os.path.join(output_dir, f"{datestamp_min}_geneformer_cellClassifier_{args.output_prefix}"))
            print(f"copying the finetuned model from {os.path.join(args.tmp_dir, datestamp, f'{datestamp_min}_geneformer_cellClassifier_{args.output_prefix}')} to {os.path.join(output_dir, f'{datestamp_min}_geneformer_cellClassifier_{args.output_prefix}')}", flush = True)
            shutil.copytree(os.path.join(args.tmp_dir, datestamp, f"{datestamp_min}_geneformer_cellClassifier_{args.output_prefix}"), os.path.join(output_dir, f"{datestamp_min}_geneformer_cellClassifier_{args.output_prefix}"))

    ######################
    # evaluate the model #
    ######################
    eval_skipped = False # a flag detecting if evaluation was skipped due to existing evaluation results
    # determine eval dataset
    if args.special_eval:
        eval_data_name = os.path.basename(args.eval_data_path) # define the evaluation data path
        eval_data_parent_dir = os.path.dirname(args.eval_data_path)
        special_eval_suffix = "special_eval_" # used to differentiate between special and regular evaluation in output files
    else:
        eval_data_name = f"{args.output_prefix}_labeled_test.dataset"
        eval_data_parent_dir = output_dir 
        special_eval_suffix = ""

    if not args.skip_evaluation: # skip the evaluation if --skip_evaluation is set
        # check if evaluation exists
        eval_result_path = f"{output_dir}/{args.output_prefix}_{special_eval_suffix}test_metrics_dict.pkl"
        if os.path.exists(eval_result_path) and os.path.getsize(eval_result_path) > 0 and not args.force_eval:
            print(f"Evaluation results exist: {eval_result_path}", flush = True)
            print(f"...Evaluation file size: {os.path.getsize(eval_result_path) / (1024 * 1024):.2f} MB", flush=True)
            print("...skipping evaluation and jumping to plotting confusion matrix", flush = True)
            eval_skipped = True
        else:
            print("evaluating model", flush = True)
            if args.use_tmp: # if use_tmp is set, use the /tmp directory to store a copy of the input data for faster I/O
                # copy the input data to the /tmp directory
                print(f"using {args.tmp_dir} directory to store input/test data", flush = True)
                tmp_dir = os.path.join(args.tmp_dir, datestamp) # store the input data in a temporary directory with the datestamp
                if os.path.exists(os.path.join(tmp_dir, eval_data_name)):
                    shutil.rmtree(os.path.join(tmp_dir, eval_data_name))
                print(f"copying {os.path.join(eval_data_parent_dir, eval_data_name)} to {os.path.join(tmp_dir, eval_data_name)}", flush = True)
                shutil.copytree(f"{eval_data_parent_dir}/{eval_data_name}", os.path.join(tmp_dir, eval_data_name))
                if os.path.exists(os.path.join(tmp_dir, f"{args.output_prefix}_id_class_dict.pkl")):
                    os.remove(os.path.join(tmp_dir, f"{args.output_prefix}_id_class_dict.pkl"))
                print(f"copying {os.path.join(output_dir, f'{args.output_prefix}_id_class_dict.pkl')} to {os.path.join(tmp_dir, f'{args.output_prefix}_id_class_dict.pkl')}", flush = True)
                shutil.copy(f"{output_dir}/{args.output_prefix}_id_class_dict.pkl", os.path.join(tmp_dir, f"{args.output_prefix}_id_class_dict.pkl"))
                output_dir = tmp_dir

            # prepare data for special evaluation (createing a feature named "label" in the arrow dataset)
            if args.special_eval:
                print("Preparing data for special evaluation", flush = True)
                start_time = time.time()
                import pickle
                from datasets import load_from_disk
                # load the id_class_dict and invert it
                with open(os.path.join(tmp_dir, f'{args.output_prefix}_id_class_dict.pkl'), 'rb') as f:
                    id_class_dict = pickle.load(f)
                value_to_id = {v: k for k, v in id_class_dict.items()}
                # load the dataset
                dataset = load_from_disk(args.eval_data_path)
                # rename the column
                dataset = dataset.rename_column(json.loads(args.cell_state_dict)["state_key"], "label")
                # Convert the label column to string to avoid type mismatches
                dataset = dataset.map(lambda x: {"label": str(x["label"])})
                # map the label to the class id and remove invalid translations
                def translate_to_class_id(example):
                    value = example['label']
                    class_id = value_to_id.get(value, -1)  # use -1 as an indicator for missing mapping
                    example['label'] = class_id
                    return example
                dataset = dataset.map(translate_to_class_id)
                dataset = dataset.filter(lambda x: x['label'] != -1)
                # print out unique labels left
                #print(value_to_id)
                #print(set(dataset.unique("label")))

                # save the dataset to disk using a temporary directory name
                temp_dataset_path = os.path.join(output_dir, eval_data_name + "_temp")
                dataset.save_to_disk(temp_dataset_path)
                # remove the existing dataset in the output directory if it exists
                old_dataset_path = os.path.join(output_dir, eval_data_name)
                if os.path.exists(old_dataset_path):
                    shutil.rmtree(old_dataset_path)
                # rename the temporary folder to the original name
                os.rename(temp_dataset_path, old_dataset_path)

                end_time = time.time()
                print(f"Time taken to prepare data for special evaluation: {end_time - start_time} seconds", flush = True)


            # print out the dimensions of the training data and the number of classes
            _data = load_from_disk(f"{output_dir}/{eval_data_name}")
            print(f"Evaluation data rows: {_data.num_rows}, columns: {len(_data.column_names)}", flush=True)
            print(f"Number of classes: {len(_data.unique('label'))}", flush=True)
            del _data
            gc.collect()

            start_time = time.time()
            cc = Classifier(classifier="cell",
                            cell_state_dict = json.loads(args.cell_state_dict) if args.cell_state_dict else None,
                            forward_batch_size=200,
                            training_args=training_args,
                            nproc=args.n_procs,
                            ngpu=args.n_gpu)

            all_metrics_test = cc.evaluate_saved_model(
                    model_directory=model_dir_path,
                    id_class_dict_file=f"{output_dir}/{args.output_prefix}_id_class_dict.pkl",
                    test_data_file=f"{output_dir}/{eval_data_name}",
                    output_directory=output_dir,
                    output_prefix=args.output_prefix,
                )

            end_time = time.time()
            print(f"Time taken to evaluate model: {end_time - start_time} seconds", flush = True)

            # plot confusion matrix
            if not args.skip_confusion_plot:
                print("plotting confusion matrix", flush = True)
                cc.plot_conf_mat(
                        conf_mat_dict={"Geneformer": all_metrics_test["conf_matrix"]},
                        output_directory=output_dir,
                    output_prefix=args.output_prefix,
                    custom_class_order=args.confu_plot_class_order.split(',') if args.confu_plot_class_order else None,
                    height=args.confu_plot_height,
                    width=args.confu_plot_width,
                )
                # copy the confusion matrix to the output directory
                # remove the existing file if it exists
                if os.path.exists(os.path.join(output_dir, f"{args.output_prefix}_conf_matrix.pdf")):
                    os.remove(os.path.join(output_dir, f"{args.output_prefix}_conf_matrix.pdf"))
                print(f"copying {os.path.join(output_dir, f'{args.output_prefix}_conf_mat.pdf')} to {os.path.join(working_dir, f'{args.output_prefix}/{datestamp}/', f'{args.output_prefix}_conf_matrix.pdf')}", flush = True)
                shutil.copy(os.path.join(output_dir, f"{args.output_prefix}_conf_mat.pdf"), os.path.join(working_dir, f"{args.output_prefix}/{datestamp}/", f"{args.output_prefix}_conf_matrix.pdf"))
        
        # copy the evaluation results to the output directory
        if args.use_tmp and not eval_skipped:
            output_dir = os.path.join(working_dir, f"{args.output_prefix}/{datestamp}/") # reset the output directory to the original output directory

            # delete the target file if it exists
            if os.path.exists(os.path.join(output_dir, f"{args.output_prefix}_{special_eval_suffix}test_metrics_dict.pkl")):
                os.remove(os.path.join(output_dir, f"{args.output_prefix}_{special_eval_suffix}test_metrics_dict.pkl"))

            # copy the test metrics dict to the output directory
            print(f"copying {os.path.join(args.tmp_dir, datestamp, f'{args.output_prefix}_test_metrics_dict.pkl')} to {os.path.join(output_dir, f'{args.output_prefix}_{special_eval_suffix}test_metrics_dict.pkl')}", flush = True)
            shutil.copy(os.path.join(args.tmp_dir, datestamp, f'{args.output_prefix}_test_metrics_dict.pkl'), os.path.join(output_dir, f'{args.output_prefix}_{special_eval_suffix}test_metrics_dict.pkl'))

            # delete the target file if it exists
            if os.path.exists(os.path.join(output_dir, f"{args.output_prefix}_{special_eval_suffix}pred_dict.pkl")):
                os.remove(os.path.join(output_dir, f"{args.output_prefix}_{special_eval_suffix}pred_dict.pkl"))
            # copy the pred dict to the output directory
            print(f"copying {os.path.join(args.tmp_dir, datestamp, f'{args.output_prefix}_pred_dict.pkl')} to {os.path.join(output_dir, f'{args.output_prefix}_{special_eval_suffix}pred_dict.pkl')}", flush = True)
            shutil.copy(os.path.join(args.tmp_dir, datestamp, f'{args.output_prefix}_pred_dict.pkl'), os.path.join(output_dir, f'{args.output_prefix}_{special_eval_suffix}pred_dict.pkl'))

    # delete the temporary directory if use_tmp is set
    if args.use_tmp and os.path.exists(os.path.join(args.tmp_dir, datestamp)):
        shutil.rmtree(os.path.join(args.tmp_dir, datestamp))