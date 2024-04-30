#!/bin/bash
#SBATCH -p sched_mem1TB_centos7
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH -o mainout.txt
#SBATCH -e mainerr.txt
#SBATCH --mem=4000

bash snakemakeslurm.sh
echo Done!!!
