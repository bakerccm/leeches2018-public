Running occupancy models

(1) Data file Ailaoshan_OTU_table.Rdata should already be updated.

(2) Run Ailaoshan_occupancy_MCMC_prepare.Rmd.

This generates Ailaoshan_occupancy_MCMC_prepare.Rdata and also the data files for individual models, i.e. the Ailaoshan_occupancy_MCMC_model???_?SU_data.Rdata.

(3) The .jags files should be prepared alongside Ailaoshan_occupancy_MCMC_prepare.Rmd, since the JAGS code references data defined in that code.

(4) To run the models on the Cannon cluster, the following code is required in addition to the .Rdata files and .jags files:

Ailaoshan_occupancy_MCMC_model_parallel_primary.sh
Ailaoshan_occupancy_MCMC_model_parallel_secondary.sh
Ailaoshan_occupancy_MCMC_model_parallel.R

OR (NOT USED -- SHOULD WORK BUT MAY REQUIRE ADJUSTMENT OF TIME LIMITS ETC.):
Ailaoshan_occupancy_MCMC_model_serial_primary.sh
Ailaoshan_occupancy_MCMC_model_serial_secondary.sh
Ailaoshan_occupancy_MCMC_model_serial.R

AND:
Ailaoshan_occupancy_MCMC_model_postprocess_dictab_primary.sh
Ailaoshan_occupancy_MCMC_model_postprocess_dictab_secondary.sh
Ailaoshan_occupancy_MCMC_model_postprocess_dictab_secondary.R
Ailaoshan_occupancy_MCMC_model_postprocess_dictab_summary.sh
Ailaoshan_occupancy_MCMC_model_postprocess_dictab_summary.R

(5) move these files to the cluster as follows:

	cd /Users/chris/Documents/work/leeches/leeches2018/analysis/Ailaoshan_occupancy_MCMC

	MYDEST="cbaker@login.rc.fas.harvard.edu:/n/piercefs/protected/Users/cbaker/leeches"

	scp Ailaoshan_occupancy_MCMC_model*.jags $MYDEST

	scp Ailaoshan_occupancy_MCMC_model???_?SU_data.Rdata $MYDEST

	scp Ailaoshan_occupancy_MCMC_model_parallel* $MYDEST

	scp Ailaoshan_occupancy_MCMC_model_postprocess_dictab* $MYDEST

	unset MYDEST

(6) log on to cluster and run scripts:

	cd /n/piercefs/protected/Users/cbaker/leeches

	sbatch Ailaoshan_occupancy_MCMC_model_parallel_primary.sh

This will run all the scripts listed in (4). Total wall time should be ~12 h.

(7) After verifying that code has executed correctly, transfer DIC results back from cluster:

	MYDEST="/Users/chris/Documents/work/leeches/leeches2018/analysis/Ailaoshan_occupancy_MCMC/"

	MYSOURCE="cbaker@login.rc.fas.harvard.edu:/n/piercefs/protected/Users/cbaker/leeches/Ailaoshan_occupancy_MCMC_model_postprocess_dictab_summary.txt"

	scp $MYSOURCE $MYDEST

	unset MYSOURCE MYDEST

(8) To post-process modelling results, use an interactive job on the cluster. But first, transfer the pre-MCMC data file so that this can be loaded in R alongside the modelling results:

	MYSOURCE="/Users/chris/Documents/work/leeches/leeches2018/analysis/Ailaoshan_occupancy_MCMC/Ailaoshan_occupancy_MCMC_prepare.Rdata"

	MYDEST="cbaker@login.rc.fas.harvard.edu:/n/piercefs/protected/Users/cbaker/leeches/"

	scp $MYSOURCE $MYDEST

	unset MYSOURCE MYDEST

Now use the code in Ailaoshan_occupancy_postprocess_cluster.R to process the model results in R.

See preamble in Ailaoshan_occupancy_MCMC_model_postprocess_cluster.R for details of how to launch a suitable interactive job and execute commands in R.

(9) Transfer post-processed model results back from cluster:

	MYSOURCE="cbaker@login.rc.fas.harvard.edu:/n/piercefs/protected/Users/cbaker/leeches"

	MYDEST="/Users/chris/Documents/work/leeches/leeches2018/analysis/"

	scp "$MYSOURCE/Ailaoshan_occupancy_postprocess_cluster.Rdata" "$MYDEST/"

	scp "$MYSOURCE/FigS3*_env_covariates_model*.pdf" "$MYDEST/figures"

	scp "$MYSOURCE/Ailaoshan_occupancy_MCMC_Rhat_*.pdf" "$MYDEST/Ailaoshan_occupancy_MCMC"

	unset MYSOURCE MYDEST
	
[Are there other files, like the model output .Rds files, that we should retain? These could be useful but there's about 1.5G data for each model variant.]

(10) Update Figure S3 locally.

Run more code locally ...
