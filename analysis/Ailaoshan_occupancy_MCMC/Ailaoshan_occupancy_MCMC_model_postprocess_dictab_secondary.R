# grabs DIC results from a single model

library("tidyverse")
library("AICcmodavg")
library("R2jags")

rm(list=ls())

setwd("/n/piercefs/protected/Users/cbaker/leeches/") ####

# get file names from command line arguments
args <- commandArgs(trailingOnly = TRUE)
results.filename <- args[1]    # loads output from model
output.filename <- args[2]  # file to save model output to as data stream

model.output <- readRDS(results.filename)

class(model.output) <- "rjags"

results <- tibble(
	model = NA,
	pD = DIC(model.output, return.pD = TRUE),
	DIC = DIC(model.output, return.pD = FALSE),
	Delta_DIC = NA,
	ModelLik = NA,
	DICWt = NA,
	Cum.Wt = NA
)

saveRDS(results, file = output.filename)
