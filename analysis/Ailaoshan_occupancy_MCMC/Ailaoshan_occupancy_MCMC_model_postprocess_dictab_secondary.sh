#!/bin/bash
#SBATCH -n 1  # number of cores
#SBATCH -N 1  # ensure that all cores are on one machine
#SBATCH -t 0-00:05  # runtime in D-HH:MM
#SBATCH -p pierce  # partition to submit to
#SBATCH --mem=8G  # memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=FAIL  # type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=bakerccm@gmail.com  # email to which notifications will be sent

module load R/3.5.1-fasrc01
export R_LIBS_USER=~/apps/R:$R_LIBS_USER

# $1 should be model output .Rds filename
# $2 should be output filename

Rscript Ailaoshan_occupancy_MCMC_model_postprocess_dictab_secondary.R $1 $2
