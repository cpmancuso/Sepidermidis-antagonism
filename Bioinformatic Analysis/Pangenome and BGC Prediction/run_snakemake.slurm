#!/bin/bash
#SBATCH --job-name snakemake.main
#SBATCH -p sched_mit_tami,sched_mit_orcd
#SBATCH -n 1
#SBATCH --time=2:00:00
#SBATCH --mem=8GB
#SBATCH -o mainerr.txt
#SBATCH -e mainerr.txt
#SBATCH --mail-user=YOUR_EMAIL_HERE
#SBATCH --mail-type=ALL
#SBATCH --exclude=node327

# Activate conda environment (may need to change name of env)
source activate snakemake

snakemake --profile .

# Print "Done!!!" at end of main log file
echo Done!!!
