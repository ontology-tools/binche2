#!/bin/bash
#SBATCH --job-name=find_missing_chebis
#SBATCH --output=logs/out/find_missing_chebis%A_%a.out
#SBATCH --error=logs/err/find_missing_chebis%A_%a.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --time=2:00:00
#SBATCH --account metabolinkai

echo Running on `hostname` at `date`

cd /idiap/temp/mcarlsson/binche2

source bincheEnv/bin/activate

python wikidata/find_missing_chebis.py

echo Done!