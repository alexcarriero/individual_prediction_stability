#!/bin/bash
#
#SBATCH --job=sc_2
#SBATCH --output="./logs/cv5_iter_%A_%a.out"
#
#SBATCH --ntasks=1
#SBATCH --time=60:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir="./"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a.j.carriero@umcutrecht.nl
#
#SBATCH --array=1-100

/home/julius_bs/acarriero/R-4.6.0/bin/Rscript "./execute.R" $SLURM_ARRAY_TASK_ID "cross_validation_5"

