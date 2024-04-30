#!/bin/bash
#SBATCH -p sched_mit_tami
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH -o logs/roaryout.txt
#SBATCH -e logs/roaryerr.txt
#SBATCH --mem=16000

source activate roary
/nfs/pool002/users/mancusoc/miniconda3/envs/roary/bin/bp_genbank2gff3.pl data/sepi_*/*.gbff -o roary_input_test
roary -p 8 roary_input_test/*.gff -f roary_output_test -e --mafft