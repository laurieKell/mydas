#!/bin/bash
#
#SBATCH --array=0-0
#SBATCH --job-name=randomD
#SBATCH --output=slurm_%a.out
/usr/lib/R/bin/Rscript --vanilla slurm_run.R
