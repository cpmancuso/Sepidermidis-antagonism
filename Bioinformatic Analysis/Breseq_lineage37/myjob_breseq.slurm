#!/bin/bash
#SBATCH -p sched_mit_hill,sched_mit_tami,sched_mit_chisholm
#SBATCH -n 1
#SBATCH --time=12:00:00
#SBATCH -o mainout.txt
#SBATCH -e mainerr.txt
#SBATCH --mem=4000

bash snakemakeslurm.sh
echo Done!!!
