#!/bin/bash
#SBATCH -J dictab_primary  # job name
#SBATCH -n 1  # number of cores
#SBATCH -N 1  # ensure that all cores are on one machine
#SBATCH -t 0-00:02  # runtime in D-HH:MM
#SBATCH -p pierce # partition to submit to
#SBATCH --mem=1M  # memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=FAIL  # type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=bakerccm@gmail.com  # email to which notifications will be sent
#SBATCH -o Ailaoshan_occupancy_MCMC_model_postprocess_dictab_primary.out
#SBATCH -e Ailaoshan_occupancy_MCMC_model_postprocess_dictab_primary.err

cd /n/piercefs/protected/Users/cbaker/leeches

# process .Rds output file from each model and store DIC results in *dictab.Rds

databases=( "LSU" "SSU" )
models=( "model01a" "model01b" "model01c" "model12a" "model12b" "model12c" "model23a" "model23b" "model23c" "model18a" "model18b" "model18c" )
jobIDs=""

for d in ${databases[@]}; do
  for m in ${models[@]}; do
    # launch job and capture jobID
    # note: set -J, -o, -e here because variables cannot be passed to #SBATCH lines in secondary script
      new_job=`sbatch -J ${m}_${d}_dictab \
        -o Ailaoshan_occupancy_MCMC_${m}_${d}_dictab.out \
        -e Ailaoshan_occupancy_MCMC_${m}_${d}_dictab.err \
        Ailaoshan_occupancy_MCMC_model_postprocess_dictab_secondary.sh Ailaoshan_occupancy_MCMC_${m}_${d}_parallel_output.Rds Ailaoshan_occupancy_MCMC_${m}_${d}_dictab.Rds`
    # extract just the jobID number and append to list
      new_job=${new_job#"Submitted batch job "}
      jobIDs=$jobIDs,$new_job
    # write job number to summary file
      echo $d $m job $new_job submitted
  done
done

jobIDs="${jobIDs#,}"  # remove leading delimiter from comma-separated list of jobs

# summarize runtime info for dictab scripts once all the individual scripts have finished running

sleep 5
summary_filename="Ailaoshan_occupancy_MCMC_model_postprocess_dictab_jobsummary.txt"
sacct -j $jobIDs --format=JobID,JobName%24,CPUTime,Elapsed,MaxRSS,State --units=G >>${summary_filename}

# produce summary table for DIC results

sbatch --dependency=afterany:${jobIDs/,/:} Ailaoshan_occupancy_MCMC_model_postprocess_dictab_summary.sh

## END of MCMC scripts
