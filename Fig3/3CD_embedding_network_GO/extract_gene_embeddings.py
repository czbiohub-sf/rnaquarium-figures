from geneformer import EmbExtractor
import os, pickle, argparse
import sys, logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
    force=True,  # <- ensure our config wins
)
logger = logging.getLogger("emb-extract")

args = argparse.ArgumentParser()
args.add_argument("--dry_run", action="store_true")
args.add_argument("--layer", type=int, default=-1, help="layer to extract embeddings from, -1: second to last layer, 0: last layer")
args.add_argument("--model_path_RQ", type=str, help="path to RQ model")
args.add_argument("--model_path_random", type=str, help="path to random model")
args.add_argument("--model_path_scrambled", type=str, help="path to scrambled model")
args.add_argument("--model_path_zebrahub", type=str, help="path to zebrahub model")
args.add_argument("--model_path_zebrahub_finetune", type=str, help="path to zebrahub finetune model")
args = args.parse_args()

# layer to extract embeddings from
assert args.layer in [-1, 0], "layer must be -1 or 0"
layer_name = "last_layer" if args.layer == 0 else "second_to_last_layer"

Geneformer_dir = "/path/to/Geneformer_RQfork"

# model list
model_list = [
    {
        "model_name": "RNAquarium",
        "token_dict_path": os.path.join(Geneformer_dir, "tokenize/RNAquarium/output/token_dictionaries/token_dictionary_RNAquarium.pickle"),
        "model_path": args.model_path_RQ,
        "dataset_path": os.path.join(Geneformer_dir, "tokenize/RNAquarium/output/tokenized_data/RNAquarium_bulk.dataset")
    },
    # {
    #     "model_name": "RNAquarium_initial_weights",
    #     "token_dict_path": os.path.join(Geneformer_dir, "tokenize/RNAquarium/output/token_dictionaries/token_dictionary_RNAquarium.pickle"),
    #     "model_path": os.path.join(Geneformer_dir, "pretrain/initial_weights/models/250301_092214_geneformer_RNAquarium_initial_weights_L6_emb256_SL2048_E30_B12_LR0.001_LSlinear_WU10000_Oadamw_DS1/models"),
    #     "dataset_path": os.path.join(Geneformer_dir, "tokenize/RNAquarium/output/tokenized_data/RNAquarium_bulk.dataset")
    # },
    {
        "model_name": "RNAquarium_scrambled",
        "token_dict_path": os.path.join(Geneformer_dir, "tokenize/RNAquarium/output/token_dictionaries/token_dictionary_RNAquarium.pickle"),
        "model_path": args.model_path_scrambled,
        "dataset_path": os.path.join(Geneformer_dir, "tokenize/RNAquarium_scrambled/output/tokenized_data/RNAquarium_bulk_scrambled.dataset")
    },
    {
        "model_name": "RNAquarium_randomized",
        "token_dict_path": os.path.join(Geneformer_dir, "tokenize/RNAquarium/output/token_dictionaries/token_dictionary_RNAquarium.pickle"),
        "model_path": args.model_path_random,
        "dataset_path": os.path.join(Geneformer_dir, "tokenize/RNAquarium_randomized/output/tokenized_data/RNAquarium_bulk_randomized.dataset")
    },
    {
        "model_name": "RNAquarium_zebrahubData",
        "token_dict_path": os.path.join(Geneformer_dir, "tokenize/RNAquarium/output/token_dictionaries/token_dictionary_RNAquarium.pickle"),
        "model_path": args.model_path_zebrahub,
        "dataset_path": os.path.join(Geneformer_dir, "tokenize/zebrahub_100-3000/output/tokenized_data/zebrahub_devtissue_100-3000.dataset")
    },
    {
        "model_name": "RNAquarium_zebrahub_finetune",
        "token_dict_path": os.path.join(Geneformer_dir, "tokenize/RNAquarium/output/token_dictionaries/token_dictionary_RNAquarium.pickle"),
        "model_path": args.model_path_zebrahub_finetune,
        "dataset_path": os.path.join(Geneformer_dir, "/path/to/Geneformer_RQfork/finetune/input2048_30epochs_zebrahub_100-3000_hyperopted/scan_devtissue/zebrahub_devtissue_classifier_freeze6layers_testsize20perc/25Aug14t10pm1/zebrahub_devtissue_classifier_freeze6layers_testsize20perc_labeled_train.dataset")
    }
]

# output directory
output_dir = "./extracted_embeddings"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# extract embeddings for each item in item_list
for model in model_list:
    token_dict_path = model["token_dict_path"]
    model_path = model["model_path"]
    model_name = model["model_name"]
    dataset_path = model["dataset_path"]
    # load token dict
    with open(token_dict_path, 'rb') as file:
        token_dict = pickle.load(file)
            
    print(f"extracting embeddings for {model_name}")
    # initiate EmbExtractor
    embex = EmbExtractor(model_type="Pretrained", # options: "Pretrained", "GeneClassifier", "CellClassifier"
                        emb_layer=args.layer, # default = -1 (second to last layer), 0: last layer
                        emb_mode="gene",
                        #emb_label=["disease","cell_type"],
                        #labels_to_plot=["disease"],
                        forward_batch_size=200,
                        nproc=8,
                        token_dictionary_file=token_dict_path) # change from current default dictionary for 30M model series

    # extracts embedding from input data
    if not args.dry_run:
        embs = embex.extract_embs(model_path, # example 30M fine-tuned model
                                dataset_path,
                                output_dir,
                                f"{model_name}_{layer_name}")
        print(f"{layer_name} embeddings extracted for {model_name}")
    else:
        logger.info(f"Dry run for {model_name}")
        # check if all the paths exist
        logger.info(f"...Checking if token_dict_path, model_path, and dataset_path exist...")
        flag = False
        if not os.path.exists(token_dict_path):
            logger.error(f"token_dict_path does not exist: {token_dict_path}")
            flag = True
        if not os.path.exists(model_path):
            logger.error(f"model_path does not exist: {model_path}")
            flag = True
        if not os.path.exists(dataset_path):
            logger.error(f"dataset_path does not exist: {dataset_path}")
            flag = True
        if flag:
            print(f"...file check failed")
        else:
            print(f"...ok")


