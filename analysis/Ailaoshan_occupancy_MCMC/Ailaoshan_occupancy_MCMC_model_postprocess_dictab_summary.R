# collates all dic results

library("tidyverse")
library("AICcmodavg")
library("R2jags")

rm(list=ls())

setwd("/n/piercefs/protected/Users/cbaker/leeches/") ####

my.databases <- c("LSU", "SSU")
my.models <- c("model01a", "model01b", "model01c", "model12a", "model12b", "model12c",  "model23a", "model23b", "model23c", "model18a", "model18b", "model18c")

model.output <- list()
for (d in my.databases) {
	model.output[[d]] <- list()
	for (m in my.models) {
		model.output[[d]][[m]] <- readRDS(paste0("Ailaoshan_occupancy_MCMC_", m , "_", d, "_dictab.Rds"))
		model.output[[d]][[m]]$model <- m
	}
}

my.dictab <- list(
	LSU = do.call(rbind, model.output[["LSU"]]) %>% as_tibble(),
	SSU = do.call(rbind, model.output[["SSU"]]) %>% as_tibble()
)

for (d in my.databases) {
	my.dictab[[d]]$Delta_DIC <- my.dictab[[d]]$DIC - min(my.dictab[[d]]$DIC)
    my.dictab[[d]]$ModelLik <- exp(-0.5 * my.dictab[[d]]$Delta_DIC)
    my.dictab[[d]]$DICWt <- my.dictab[[d]]$ModelLik/sum(my.dictab[[d]]$ModelLik)
    my.dictab[[d]] <- my.dictab[[d]][order(my.dictab[[d]]$Delta_DIC), ]
    my.dictab[[d]]$Cum.Wt <- cumsum(my.dictab[[d]]$DICWt)
}

cat("\nLSU\n")
my.dictab[["LSU"]]
cat("\nSSU\n")
my.dictab[["SSU"]]
