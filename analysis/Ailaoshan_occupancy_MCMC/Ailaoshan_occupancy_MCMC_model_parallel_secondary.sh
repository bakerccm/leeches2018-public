#!/bin/bash
#SBATCH -n 3  # number of cores
#SBATCH -N 1  # ensure that all cores are on one machine
#SBATCH -t 0-08:00  # runtime in D-HH:MM
#SBATCH -p pierce  # partition to submit to
#SBATCH --mem=24G  # memory pool for all cores (see also --mem-per-cpu)

# note: -J, -o and -e set in primary .sh file so variables can be passed to them

# $1 should be modelxx
# $2 should be dataset (i.e. LSU or SSU)

module load R/3.5.1-fasrc01
export R_LIBS_USER=~/apps/R:$R_LIBS_USER

Rscript Ailaoshan_occupancy_MCMC_model_parallel.R Ailaoshan_occupancy_MCMC_$1_$2_data.Rdata Ailaoshan_occupancy_MCMC_$1.jags Ailaoshan_occupancy_MCMC_$1_$2_parallel_output.Rds

sleep 5
sacct -j $SLURM_JOBID --format=JobID,JobName%24,CPUTime,Elapsed,MaxRSS --units=G
