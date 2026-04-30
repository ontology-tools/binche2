#!/bin/bash
#SBATCH --job-name=extract_hmdb
#SBATCH --output=logs/out/extract_hmdb%A_%a.out
#SBATCH --error=logs/err/extract_hmdb%A_%a.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --time=5:00:00
#SBATCH --account metabolinkai

echo Running on `hostname` at `date`

cd /idiap/temp/mcarlsson/binche2

source bincheEnv/bin/activate

python hmdb/extract_hmdb.py

echo Done!