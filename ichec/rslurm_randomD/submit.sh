#!/bin/bash
#
#SBATCH --array=0-0
#SBATCH --job-name=randomD
#SBATCH --output=slurm_%a.out

#SBATCH --nodes=1

#SBATCH --time=00:20:00
/usr/lib/R/bin/Rscript --vanilla slurm_run.R
