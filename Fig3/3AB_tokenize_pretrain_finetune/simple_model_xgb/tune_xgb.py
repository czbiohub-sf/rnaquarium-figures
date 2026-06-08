import os
import pickle
import argparse
import tempfile
import numpy as np
import xgboost as xgb
from typing import Dict, Any
import matplotlib.pyplot as plt
import ray
from ray import tune
from ray.tune.schedulers import ASHAScheduler
from ray.tune.search.optuna import OptunaSearch
from sklearn.metrics import auc, roc_curve, f1_score  

def load_data(file_path: str) -> Dict[str, np.ndarray]:
    """
    Load data from pickle file
    
    Args:
        file_path: Path to the pickle file
        
    Returns:
        Dictionary containing 'X' and 'y' numpy arrays
    """
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
    return data


def train_xgboost(config: Dict[str, Any], train_data: Dict[str, np.ndarray], 
                  eval_data: Dict[str, np.ndarray], objective: str, 
                  eval_metrics: list, use_gpu: bool = False):
    """
    Training function for XGBoost model that Ray Tune will call
    """
    # Convert data to DMatrix format
    dtrain = xgb.DMatrix(train_data['X'], label=train_data['y'])
    deval = xgb.DMatrix(eval_data['X'], label=eval_data['y'])
    
    # Set up parameters
    param = {
        'objective': objective,
        'eval_metric': eval_metrics,
        'max_depth': int(config['max_depth']),
        'learning_rate': config['learning_rate'],
        'subsample': config['subsample'],
        'colsample_bytree': config['colsample_bytree'],
        'min_child_weight': config['min_child_weight'],
        'alpha': config['alpha'],
        'lambda': config['lambda'],
        'gamma': config['gamma'],
        'seed': 42,
    }
    
    # Add GPU configuration if needed
    if use_gpu:
        param['device'] = 'cuda'
    else:
        param['device'] = 'cpu'

    # Automatically add num_class for multi-class objectives
    if objective.startswith("multi"):
        num_class = int(np.unique(train_data['y']).shape[0])
        param['num_class'] = num_class

    # Train model
    evals_result = {}
    bst = xgb.train(
        param,
        dtrain,
        num_boost_round=int(config['num_boost_round']),
        evals=[(dtrain, 'train'), (deval, 'eval')],
        evals_result=evals_result,
        early_stopping_rounds=50,
        verbose_eval=False
    )
    
    # Report metrics to Ray Tune
    tune_report_dict = {}
    for m in eval_metrics:
        tune_report_dict[f"eval_{m}"] = evals_result['eval'][m][-1]
    tune_report_dict["best_iteration"] = bst.best_iteration
    
    tune.report(tune_report_dict)


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='XGBoost hyperparameter tuning with Ray')
    
    # Data parameters
    parser.add_argument('--train_data', type=str, default='data/OrgIP_train_data.pkl', 
                        help='Path to training data pickle file (default: OrgIP_train_data.pkl)')
    parser.add_argument('--eval_data', type=str, default='data/OrgIP_eval_data.pkl', 
                        help='Path to evaluation data pickle file (default: OrgIP_eval_data.pkl)')
    
    # XGBoost parameters
    parser.add_argument('--objective', type=str, default='multi:softprob', 
                        help='XGBoost objective function (default: multi:softprob). '
                             'Options: binary:logistic, reg:squarederror, multi:softmax, etc.')
    parser.add_argument('--eval_metrics', type=str, default='logloss,error', 
                        help='Comma-separated list of evaluation metrics (default: logloss,error)')
    
    # Hardware parameters
    parser.add_argument('--use_gpu', action='store_true', 
                        help='Use GPU for training (default: False)')
    parser.add_argument('--cpus_per_trial', type=int, default=1, 
                        help='CPUs per trial (default: 1)')
    parser.add_argument('--gpus_per_trial', type=float, default=0.1, 
                        help='GPUs per trial when using GPU (default: 0.1)')
    
    # Tuning parameters
    parser.add_argument('--num_trials', type=int, default=100, 
                        help='Number of hyperparameter configurations to try (default: 100)')
    parser.add_argument('--max_num_epochs', type=int, default=100, 
                        help='Maximum number of training epochs (default: 100)')
    parser.add_argument('--concurrent_trials', type=int, default=10,
                        help='Number of trials to run concurrently (default: 10)')

    # Output parameters
    parser.add_argument('--output_dir', type=str, default='./output', 
                        help='Directory to save output files (default: ./output)')
    
    args = parser.parse_args()


     # Create output directory for final results
    os.makedirs(args.output_dir, exist_ok=True)
    # Create temporary directory for Ray Tune to use
    temp_dir = tempfile.mkdtemp()
    print(f"Using temporary directory for Ray Tune: {temp_dir}")
    
    ############
    # Load data
    print(f"Loading training data from {args.train_data}")
    train_data = load_data(args.train_data)
    print(f"Loading evaluation data from {args.eval_data}")
    eval_data = load_data(args.eval_data)

    # Check if eval_data['y'] is a numpy array; if not, convert it.
    if not isinstance(eval_data['y'], np.ndarray):
        eval_data['y'] = np.array(eval_data['y'])
    
    # Parse evaluation metrics
    eval_metrics = args.eval_metrics.split(',')
    if args.objective.startswith("multi"):
        # Check if user hasn't provided custom metrics different from defaults
        if eval_metrics == ['logloss', 'error']:
            eval_metrics = ['mlogloss', 'merror']
    print(f"Using objective: {args.objective}")
    print(f"Using evaluation metrics: {eval_metrics}")
    
    # Initialize Ray
    print("Initializing Ray...")
    ray.init()
    
    # Define search space
    search_space = {
        'max_depth': tune.qrandint(3, 15, 1),
        'learning_rate': tune.loguniform(0.001, 0.5),
        'subsample': tune.uniform(0.5, 1.0),
        'colsample_bytree': tune.uniform(0.5, 1.0),
        'min_child_weight': tune.qrandint(1, 10, 1),
        'alpha': tune.loguniform(1e-8, 1.0),
        'lambda': tune.loguniform(1e-8, 1.0),
        'gamma': tune.loguniform(1e-8, 1.0),
        'num_boost_round': tune.qrandint(50, 500, 10)
    }

    # Define Optuna search algorithm
    search_alg = OptunaSearch()
    
    # Define scheduler
    scheduler = ASHAScheduler(
        max_t=args.max_num_epochs,
        grace_period=10,
        reduction_factor=2
    )
    
    # Resources per trial
    if args.use_gpu:
        resources_per_trial = {"gpu": args.gpus_per_trial}
        print(f"Using GPU mode with {args.gpus_per_trial} GPUs per trial")
    else:
        resources_per_trial = {"cpu": args.cpus_per_trial}
        print(f"Using CPU mode with {args.cpus_per_trial} CPUs per trial")
    
    print(f"Starting hyperparameter tuning with {args.num_trials} trials...")
    
    ###########################
    # Run hyperparameter tuning
    analysis = tune.run(
        tune.with_parameters(train_xgboost, 
                            train_data=train_data, 
                            eval_data=eval_data, 
                            objective=args.objective,
                            eval_metrics=eval_metrics,
                            use_gpu=args.use_gpu),
        resources_per_trial=resources_per_trial,
        metric=f"eval_{eval_metrics[0]}",
        mode="min",  # Assuming first metric should be minimized
        config=search_space,
        num_samples=args.num_trials,
        scheduler=scheduler,
        search_alg=search_alg,
        verbose=2,
        storage_path=temp_dir, 
        name="xgboost_tuning",
    )
    
    # Get best configuration
    best_config = analysis.get_best_config(metric=f"eval_{eval_metrics[0]}", mode="min")
    print("\nBest hyperparameters:")
    for param, value in best_config.items():
        print(f"  {param}: {value}")

    # Save the best config
    config_path = os.path.join(args.output_dir, 'best_config.pkl')
    with open(config_path, 'wb') as f:
        pickle.dump(best_config, f)
    config_path_txt = os.path.join(args.output_dir, 'best_config.txt')
    with open(config_path_txt, 'w') as f:
        for param, value in best_config.items():
            f.write(f"{param}: {value}\n")
    print(f"Best configuration saved as '{config_path}' and '{config_path_txt}'")
    
    # Display best trial results
    best_trial = analysis.best_trial
    print("\nBest trial metrics:")
    for metric, value in best_trial.last_result.items():
        if metric.startswith("eval_"):
            print(f"  {metric}: {value:.6f}")
    
    print("\nHyperparameter tuning completed successfully!")
    ray.shutdown()
    
    print("\nTraining final model with best hyperparameters...")
    
    ##########################################
    # Train a final model with the best config
    dtrain = xgb.DMatrix(train_data['X'], label=train_data['y'])
    deval = xgb.DMatrix(eval_data['X'], label=eval_data['y'])
    
    # Set up parameters
    param = {
        'objective': args.objective,
        'eval_metric': eval_metrics,
        'max_depth': int(best_config['max_depth']),
        'learning_rate': best_config['learning_rate'],
        'subsample': best_config['subsample'],
        'colsample_bytree': best_config['colsample_bytree'],
        'min_child_weight': best_config['min_child_weight'],
        'alpha': best_config['alpha'],
        'lambda': best_config['lambda'],
        'gamma': best_config['gamma'],
        'seed': 42,
    }
    
    # Add GPU configuration if needed
    if args.use_gpu:
        param['device'] = 'cuda'
    else:
        param['device'] = 'cpu'
    
    # Automatically add num_class for multi-class objectives
    if args.objective.startswith("multi"):
        num_class = int(np.unique(train_data['y']).shape[0])
        param['num_class'] = num_class

    # Train final model
    final_model = xgb.train(
        param,
        dtrain,
        num_boost_round=int(best_config['num_boost_round']),
        evals=[(dtrain, 'train'), (deval, 'eval')],
        early_stopping_rounds=50,
        verbose_eval=10
    )
    
    # Save the model
    model_path = os.path.join(args.output_dir, 'best_xgboost_model.json')
    final_model.save_model(model_path)
    print(f"Final model saved as '{model_path}'")
    
    ####################################
    # Calculate and print macro F1 score
    # Predict probabilities for each class using the final model
    y_score = final_model.predict(xgb.DMatrix(eval_data['X']))
    y_pred = np.argmax(y_score, axis=1)

    macro_f1 = f1_score(eval_data['y'], y_pred, average='macro')
    print(f"\nMacro F1 Score on evaluation data: {macro_f1:.6f}")
    metrics_dict = {'macro_f1': macro_f1}
    
    # Save the macro F1 score
    metrics_path = os.path.join(args.output_dir, 'xgb_test_metrics_dict.pkl')
    with open(metrics_path, 'wb') as f:
        pickle.dump(metrics_dict, f)

    print(f"Evaluation metrics saved as '{metrics_path}'")

    #######################
    # Draw an AUCROC curve 

    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    for i in range(num_class):
        fpr[i], tpr[i], _ = roc_curve(eval_data['y'] == i, y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    # Define label_mapping: use train_data["int_to_label_dict"] if it exists; otherwise default to string of class index
    if train_data.get("int_to_label_dict"):
        # Invert the dictionary so that index -> label name
        label_mapping =  train_data["int_to_label_dict"]
    else:
        label_mapping = {i: str(i) for i in range(num_class)}

    # Plot the ROC curve for each class
    plt.figure(figsize=(10, 9))
    colors = plt.cm.tab20(range(20))[0:num_class]

    for i, color in zip(range(num_class), colors):
        plt.plot(fpr[i], tpr[i], color=color, lw=2,
                label='class {0} ({1}), area = {2:0.4f}'
                .format(i, label_mapping[i], roc_auc[i])
                )

    plt.plot([0, 1], [0, 1], 'k--', lw=1)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic for multi-class data')
    plt.legend(loc="lower right")
    plt.savefig(os.path.join(args.output_dir, "AUC_ROC.pdf"), format="pdf", bbox_inches='tight')
    print(f"AUC ROC curve saved as '{os.path.join(args.output_dir, 'AUC_ROC.pdf')}'")


if __name__ == "__main__":
    main()