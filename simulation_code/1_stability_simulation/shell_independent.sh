#!/bin/bash
#
#SBATCH --job=sc_1
#SBATCH --output="./logs/independent_iter_%A_%a.out"
#
#SBATCH --ntasks=1
#SBATCH --time=25:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir="./"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a.j.carriero@umcutrecht.nl
#
#SBATCH --array=1-100

/home/julius_bs/acarriero/R-4.6.0/bin/Rscript "./execute.R" $SLURM_ARRAY_TASK_ID "independent"

