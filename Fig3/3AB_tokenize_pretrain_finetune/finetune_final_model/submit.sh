# prep data
sbatch cmd_devtissue_prepdata.sh 25Aug14t10pm1
sbatch cmd_devtissue_prepdata.sh 25Aug14t10pm2
sbatch cmd_devtissue_prepdata.sh 25Aug14t10pm3
sbatch cmd_devtissue_prepdata.sh 25Aug14t10pm4
sbatch cmd_devtissue_prepdata.sh 25Aug14t10pm5
sbatch cmd_devtissue_prepdata.sh 25Aug14t10pm6


sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm1 0.05 0.05 0.05
sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm1 0.10 0.05 0.10
sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm1 0.15 0.05 0.20
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm1 0.25 0.05 0.30
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm1 0.35 0.05 0.50
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm1 0.55 0.05 0.70
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm1 0.75 0.05 0.95

sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm2 0.05 0.05 0.05
sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm2 0.10 0.05 0.10
sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm2 0.15 0.05 0.20
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm2 0.25 0.05 0.30
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm2 0.35 0.05 0.50
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm2 0.55 0.05 0.70
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm2 0.75 0.05 0.95

sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm3 0.05 0.05 0.05
sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm3 0.10 0.05 0.10
sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm3 0.15 0.05 0.20
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm3 0.25 0.05 0.30
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm3 0.35 0.05 0.50
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm3 0.55 0.05 0.70
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm3 0.75 0.05 0.95

sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm4 0.05 0.05 0.05
sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm4 0.10 0.05 0.10
sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm4 0.15 0.05 0.20
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm4 0.25 0.05 0.30
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm4 0.35 0.05 0.50
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm4 0.55 0.05 0.70
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm4 0.75 0.05 0.95

sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm5 0.05 0.05 0.05
sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm5 0.10 0.05 0.10
sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm5 0.15 0.05 0.20
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm5 0.25 0.05 0.30
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm5 0.35 0.05 0.50
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm5 0.55 0.05 0.70
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm5 0.75 0.05 0.95

sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm6 0.05 0.05 0.05
sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm6 0.10 0.05 0.10
sbatch --gpus-per-task=4 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm6 0.15 0.05 0.20
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm6 0.25 0.05 0.30
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm6 0.35 0.05 0.50
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm6 0.55 0.05 0.70
sbatch --gpus-per-task=2 cmd_devtissue_finetuneOnly.sh 25Aug14t10pm6 0.75 0.05 0.95

