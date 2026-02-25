#!/bin/bash
#SBATCH --job-name=create_files
#SBATCH --output=logs/out/create_files%A_%a.out
#SBATCH --error=logs/err/create_files%A_%a.err
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20GB
#SBATCH --gpus rtx3090:1
#SBATCH --time=20:00:00
#SBATCH --account metabolinkai

echo Running on `hostname` at `date`

cd /idiap/temp/mcarlsson/binche2

source bincheEnv/bin/activate

python create_files.py 

echo Done!