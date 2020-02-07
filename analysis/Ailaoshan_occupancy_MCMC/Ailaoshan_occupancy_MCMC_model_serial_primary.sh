#!/bin/bash
#SBATCH -J serial_primary  # job name
#SBATCH -n 1  # number of cores
#SBATCH -N 1  # ensure that all cores are on one machine
#SBATCH -t 0-00:02  # runtime in D-HH:MM
#SBATCH -p pierce  # partition to submit to
#SBATCH --mem=1M  # memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=FAIL  # type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=bakerccm@gmail.com  # email to which notifications will be sent
#SBATCH -o Ailaoshan_occupancy_MCMC_model_serial_primary.out
#SBATCH -e Ailaoshan_occupancy_MCMC_model_serial_primary.err

cd /n/piercefs/protected/Users/cbaker/leeches

summary_filename="Ailaoshan_occupancy_MCMC_model_serial_summary.txt"
echo "Start time: " `date` >${summary_filename}

# launch scripts to run individual models

databases=( "LSU" "SSU" )
models=( "model01a" "model01b" "model01c" "model12a" "model12b" "model12c" "model23a" "model23b" "model23c" "model18a" "model18b" "model18c" )
jobIDs=""

for d in ${databases[@]}; do
  for m in ${models[@]}; do
    # launch job and capture jobID
    # note: set -J, -o, -e here because variables cannot be passed to #SBATCH lines in secondary script
      new_job=`sbatch -J ${m}_${d}_serial \
        -o Ailaoshan_occupancy_MCMC_${m}_${d}_serial.out -e Ailaoshan_occupancy_MCMC_${m}_${d}_serial.err \
        Ailaoshan_occupancy_MCMC_model_serial_secondary.sh ${m} ${d}`
    # extract just the jobID number and append to list
      new_job=${new_job#"Submitted batch job "}
      jobIDs=$jobIDs,$new_job
    # write job number to summary file
      echo $d $m job $new_job submitted >>${summary_filename}
  done
done

jobIDs="${jobIDs#,}"  # remove leading delimiter from comma-separated list of jobs

# summarize runtime info once individual models are finished running

sleep 5
echo '#!/bin/bash' >serial_summary.sh
echo "sacct -j $jobIDs --format=JobID,JobName%24,CPUTime,Elapsed,MaxRSS,State --units=G >>${summary_filename}" >>serial_summary.sh
echo "echo 'End time: ' \`date\` >>${summary_filename}" >>serial_summary.sh
sbatch -J serial_summary -n 1 -N 1 -t 0-00:01 -p pierce --mem=1M \
    --mail-type=END --mail-user=bakerccm@gmail.com \
    -o serial_summary.out -e serial_summary.err \
    --dependency=afterany:${jobIDs/,/:} serial_summary.sh
rm serial_summary.sh
