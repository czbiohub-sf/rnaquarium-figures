# prep data
sbatch cmd_devtissue_prepdata.sh 25Aug09t9am1

# finetune/hyperopt
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.05 0.05 0.05
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.10 0.05 0.10
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.15 0.05 0.15
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.20 0.05 0.20
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.25 0.05 0.25
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.30 0.05 0.30
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.35 0.05 0.35
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.40 0.05 0.40
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.45 0.05 0.45
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.50 0.05 0.50
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.55 0.05 0.55
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.60 0.05 0.60
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.65 0.05 0.65
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.70 0.05 0.70
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.75 0.05 0.75
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.80 0.05 0.80
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.85 0.05 0.85
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.90 0.05 0.90
sbatch cmd_devtissue_finetuneOnly.sh 25Aug09t9am1 0.95 0.05 0.95

# process hyperopt results
