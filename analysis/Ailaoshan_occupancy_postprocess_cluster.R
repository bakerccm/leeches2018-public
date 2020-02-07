# R code for postprocessing occupancy modelling results before manuscript analysis

# This could be run as a script, but it's easier to do this interactively so we can pick the
# best models based on DIC and just extract data and figures for those.

########################################################################################
## run this on cluster ##

# (need to ssh with "-Y" argument for this to work)
ssh -Y cbaker@login.rc.fas.harvard.edu
srun --mem 16G -p test -t 06:00:00 --pty --x11=first bash
module load R/3.5.1-fasrc01
export R_LIBS_USER=~/apps/R:$R_LIBS_USER
R

########################################################################################
### now, in R ... ##

library("tidyverse")
library("ggpubr")

rm(list=ls())

setwd("/n/piercefs/protected/Users/cbaker/leeches")

# need this for zst
    load("Ailaoshan_occupancy_MCMC_prepare.Rdata")

# model outputs generated on cluster:

    # Ailaoshan_occupancy_MCMC_model01a_LSU_parallel_output.Rds
    # Ailaoshan_occupancy_MCMC_model01b_LSU_parallel_output.Rds
    # Ailaoshan_occupancy_MCMC_model01c_LSU_parallel_output.Rds

    # Ailaoshan_occupancy_MCMC_model12a_LSU_parallel_output.Rds ## LSU model no.3
    # Ailaoshan_occupancy_MCMC_model12b_LSU_parallel_output.Rds
    # Ailaoshan_occupancy_MCMC_model12c_LSU_parallel_output.Rds

    # Ailaoshan_occupancy_MCMC_model18a_LSU_parallel_output.Rds ## LSU model no.2
    # Ailaoshan_occupancy_MCMC_model18b_LSU_parallel_output.Rds
    # Ailaoshan_occupancy_MCMC_model18c_LSU_parallel_output.Rds

    # Ailaoshan_occupancy_MCMC_model23a_LSU_parallel_output.Rds ## LSU model no.1 ##
    # Ailaoshan_occupancy_MCMC_model23b_LSU_parallel_output.Rds
    # Ailaoshan_occupancy_MCMC_model23c_LSU_parallel_output.Rds

    # Ailaoshan_occupancy_MCMC_model01a_SSU_parallel_output.Rds
    # Ailaoshan_occupancy_MCMC_model01b_SSU_parallel_output.Rds
    # Ailaoshan_occupancy_MCMC_model01c_SSU_parallel_output.Rds

    # Ailaoshan_occupancy_MCMC_model12a_SSU_parallel_output.Rds ## SSU model no.3
    # Ailaoshan_occupancy_MCMC_model12b_SSU_parallel_output.Rds
    # Ailaoshan_occupancy_MCMC_model12c_SSU_parallel_output.Rds

    # Ailaoshan_occupancy_MCMC_model18a_SSU_parallel_output.Rds ## SSU model no.2
    # Ailaoshan_occupancy_MCMC_model18b_SSU_parallel_output.Rds
    # Ailaoshan_occupancy_MCMC_model18c_SSU_parallel_output.Rds

    # Ailaoshan_occupancy_MCMC_model23a_SSU_parallel_output.Rds ## SSU model no.1 ##
    # Ailaoshan_occupancy_MCMC_model23b_SSU_parallel_output.Rds
    # Ailaoshan_occupancy_MCMC_model23c_SSU_parallel_output.Rds

# get filenames for top models in each dataset

    MCMC.filenames <- list(
        LSU = list(model23a = paste0("Ailaoshan_occupancy_MCMC_", "model23a", "_", "LSU", "_parallel_output.Rds")),
        SSU = list(model23a = paste0("Ailaoshan_occupancy_MCMC_", "model23a", "_", "SSU", "_parallel_output.Rds"))
    )

########################################################################################
## export model output from top models in each dataset for local analysis ##

# parse model estimates from JAGS output
    
    z.output <- list()
    psi.output <- list()
    Nsite.output <- list()
    lpsi0.output <- list()
    lp0.output <- list()
    occ.output <- list()

    for (i in names(MCMC.filenames)) {

        z.output[[i]] <- list()
        psi.output[[i]] <- list()
        Nsite.output[[i]] <- list()
        lpsi0.output[[i]] <- list()
        lp0.output[[i]] <- list()
        occ.output[[i]] <- list()
    
        for (m in names(MCMC.filenames[[i]])) {
        
            model.output <- readRDS(file = MCMC.filenames[[i]][[m]])
        
            # z varies by species and site
                z.output[[i]][[m]] <- model.output$BUGSoutput$summary[grep("z\\[", rownames(model.output$BUGSoutput$summary)),]
                rn <- rownames(z.output[[i]][[m]])
                z.output[[i]][[m]] <- as_tibble(z.output[[i]][[m]]) %>%
                    # store results more efficiently since these have to be integer
                        mutate(`2.5%`=as.integer(`2.5%`), `25%`=as.integer(`25%`), `50%`=as.integer(`50%`),
                            `75%`=as.integer(`75%`), `97.5%`=as.integer(`97.5%`), n.eff = as.integer(n.eff)) %>%
                    # parse rownames into indexes representing Polygon_ID and OTU
                        mutate(var = rn) %>%
                        mutate(var = gsub("z\\[","",var)) %>% mutate(var = gsub("\\]","",var)) %>%
                        separate(var, sep=',', into=c("Polygon_ID_","OTU_")) %>%
                    # convert indexes back to our original Polygon_ID and OTU labels
                        mutate(Polygon_ID = rownames(zst[[i]])[as.numeric(Polygon_ID_)]) %>%
                        mutate(OTU = colnames(zst[[i]])[as.numeric(OTU_)]) %>%
                        select(-Polygon_ID_, -OTU_)
                rm(rn)

            # psi varies by species and site
                psi.output[[i]][[m]] <- model.output$BUGSoutput$summary[grep("psi\\[", rownames(model.output$BUGSoutput$summary)),]
                rn <- rownames(psi.output[[i]][[m]])
                psi.output[[i]][[m]] <- as_tibble(psi.output[[i]][[m]]) %>%
                    # store results more efficiently since these have to be integer
                        mutate(n.eff = as.integer(n.eff)) %>%
                    # parse rownames into indexes representing Polygon_ID and OTU
                        mutate(var = rn) %>%
                        mutate(var = gsub("psi\\[","",var)) %>% mutate(var = gsub("\\]","",var)) %>%
                        separate(var, sep=',', into=c("Polygon_ID_","OTU_")) %>%
                    # convert indexes back to our original Polygon_ID and OTU labels
                        mutate(Polygon_ID = rownames(zst[[i]])[as.numeric(Polygon_ID_)]) %>%
                        mutate(OTU = colnames(zst[[i]])[as.numeric(OTU_)]) %>%
                        select(-Polygon_ID_, -OTU_)
                rm(rn)

            # Nsite varies by site
                Nsite.output[[i]][[m]] <- model.output$BUGSoutput$summary[grep("Nsite\\[", rownames(model.output$BUGSoutput$summary)),]
                rn <- rownames(Nsite.output[[i]][[m]])
                Nsite.output[[i]][[m]] <- as_tibble(Nsite.output[[i]][[m]]) %>%
                    # store results more efficiently since these have to be integer
                        mutate(n.eff = as.integer(n.eff)) %>%
                    # parse rownames into indexes representing Polygon_ID
                        mutate(Polygon_ID_ = rn) %>%
                        mutate(Polygon_ID_ = gsub("Nsite\\[","",Polygon_ID_)) %>% mutate(Polygon_ID_ = gsub("\\]","",Polygon_ID_)) %>%
                        mutate(Polygon_ID_ = as.numeric(Polygon_ID_)) %>%
                    # convert indexes back to our original Polygon_ID labels
                        mutate(Polygon_ID = rownames(zst[[i]])[Polygon_ID_]) %>%
                        select(-Polygon_ID_)
                rm(rn)

            # lpsi0 varies by species
                lpsi0.output[[i]][[m]] <- model.output$BUGSoutput$summary[grep("lpsi0\\[", rownames(model.output$BUGSoutput$summary)),]
                rn <- rownames(lpsi0.output[[i]][[m]])
                lpsi0.output[[i]][[m]] <- as_tibble(lpsi0.output[[i]][[m]]) %>%
                    # store results more efficiently since these have to be integer
                        mutate(n.eff = as.integer(n.eff)) %>%
                    # parse rownames into indexes representing Polygon_ID
                        mutate(OTU_ = rn) %>%
                        mutate(OTU_ = gsub("lpsi0\\[","",OTU_)) %>% mutate(OTU_ = gsub("\\]","",OTU_)) %>%
                        mutate(OTU_ = as.numeric(OTU_)) %>%
                    # convert indexes back to our original OTU labels
                        mutate(OTU = colnames(zst[[i]])[OTU_]) %>%
                        select(-OTU_)
                rm(rn)

            # lp0 varies by species
                lp0.output[[i]][[m]] <- model.output$BUGSoutput$summary[grep("lp0\\[", rownames(model.output$BUGSoutput$summary)),]
                rn <- rownames(lp0.output[[i]][[m]])
                lp0.output[[i]][[m]] <- as_tibble(lp0.output[[i]][[m]]) %>%
                    # store results more efficiently since these have to be integer
                        mutate(n.eff = as.integer(n.eff)) %>%
                    # parse rownames into indexes representing Polygon_ID
                        mutate(OTU_ = rn) %>%
                        mutate(OTU_ = gsub("lp0\\[","",OTU_)) %>% mutate(OTU_ = gsub("\\]","",OTU_)) %>%
                        mutate(OTU_ = as.numeric(OTU_)) %>%
                    # convert indexes back to our original OTU labels
                        mutate(OTU = colnames(zst[[i]])[OTU_]) %>%
                        select(-OTU_)
                rm(rn)

            # occ varies by species
                occ.output[[i]][[m]] <- model.output$BUGSoutput$summary[grep("occ\\[", rownames(model.output$BUGSoutput$summary)),]
                rn <- rownames(occ.output[[i]][[m]])
                occ.output[[i]][[m]] <- as_tibble(occ.output[[i]][[m]]) %>%
                    # store results more efficiently since these have to be integer
                        mutate(n.eff = as.integer(n.eff)) %>%
                    # parse rownames into indexes representing Polygon_ID
                        mutate(OTU_ = rn) %>%
                        mutate(OTU_ = gsub("occ\\[","",OTU_)) %>% mutate(OTU_ = gsub("\\]","",OTU_)) %>%
                        mutate(OTU_ = as.numeric(OTU_)) %>%
                    # convert indexes back to our original OTU labels
                        mutate(OTU = colnames(zst[[i]])[OTU_]) %>%
                        select(-OTU_)
                rm(rn)
            
            rm(model.output)

        }
    }

    rm(i,m)
 
# prepare per-OTU model estimates for psi and p, converting to probability scale

    estimates.psi <- list()
    estimates.p <- list()

    for (i in names(MCMC.filenames)) {
        for (m in names(MCMC.filenames[[i]])) {
            # lpsi0 is the species-specific constant for occupancy
            # converting to probability scale should give occupancy estimate for constant (i.e. mean) values
            # for environmental covariates, since predictors are scaled and centered
            estimates.psi[[i]] <- lpsi0.output[[i]][[m]] %>%
                dplyr::select(mean, `2.5%`, `97.5%`, OTU) %>%
                rename(lpsi0_mean = mean, `lpsi0_2.5%`=`2.5%`, `lpsi0_97.5%`=`97.5%`) %>%
                mutate(prob_lpsi0_mean = boot::inv.logit(lpsi0_mean),
                       `prob_lpsi0_2.5%` = boot::inv.logit(`lpsi0_2.5%`),
                       `prob_lpsi0_97.5%` = boot::inv.logit(`lpsi0_97.5%`)) %>%
                dplyr::select(OTU, prob_lpsi0_mean, `prob_lpsi0_2.5%`, `prob_lpsi0_97.5%`)
            # lp0 is the species-specific constant for detection
            # converting to probability scale should give detection estimate for constant (i.e. mean) value
            # for numleeches, since predictors are scaled and centered
            estimates.p[[i]] <- lp0.output[[i]][[m]] %>%
                dplyr::select(mean, `2.5%`, `97.5%`, OTU) %>%
                rename(lp0_mean = mean, `lp0_2.5%`=`2.5%`, `lp0_97.5%`=`97.5%`) %>%
                mutate(prob_lp0_mean = boot::inv.logit(lp0_mean),
                       `prob_lp0_2.5%` = boot::inv.logit(`lp0_2.5%`),
                       `prob_lp0_97.5%` = boot::inv.logit(`lp0_97.5%`)) %>%
                dplyr::select(OTU, prob_lp0_mean, `prob_lp0_2.5%`, `prob_lp0_97.5%`)
        }
    }

    estimates.psi <- bind_rows(LSU = estimates.psi[["LSU"]], SSU = estimates.psi[["SSU"]], .id = "dataset")
    estimates.p <- bind_rows(LSU = estimates.p[["LSU"]], SSU = estimates.p[["SSU"]], .id = "dataset")
    OTU.estimates <- full_join(estimates.psi, estimates.p, by = c("OTU","dataset"))
    rm(estimates.psi, estimates.p)

# saves processed model.output

save(z.output, psi.output, Nsite.output, lpsi0.output, lp0.output, occ.output, OTU.estimates,
    file = "Ailaoshan_occupancy_postprocess_cluster.Rdata")

########################################################################################
## draw Fig S3 plots on cluster ##

# get data and model output

    model.data <- list()
    load("Ailaoshan_occupancy_MCMC_model23a_LSU_data.Rdata")
    model.data$LSU <- jags.data
    load("Ailaoshan_occupancy_MCMC_model23a_SSU_data.Rdata")
    model.data$SSU <- jags.data
    rm(jags.data)

    model.output <- list()
    model.output$LSU = readRDS(file = MCMC.filenames[["LSU"]][["model23a"]])
    model.output$SSU = readRDS(file = MCMC.filenames[["SSU"]][["model23a"]])

# generate predicted values

    sims <- list()
    nsamp <- list()

    elev.pred <- list(); road.pred <- list(); reserve.pred <- list(); numleeches.pred <- list()
    predC <- list()
    pmC <- list(); criI <- list()

    for (i in c("LSU","SSU")){

        # sims.list has model output arranged by variable
        # and appears to be a rearranged version of sims.matrix
        sims[[i]] <- model.output[[i]]$BUGSoutput$sims.list
        
        # number of MCMC samples
        nsamp[[i]] <- nrow(sims[[i]][[1]])
        
        # generate sequences of predictor values
        elev.pred[[i]] <- seq(min(model.data[[i]]$model.data$elev, na.rm = TRUE), max(model.data[[i]]$model.data$elev, na.rm = TRUE), length.out = 500)
        reserve.pred[[i]] <- seq(min(model.data[[i]]$model.data$reserve, na.rm = TRUE), max(model.data[[i]]$model.data$reserve, na.rm = TRUE), length.out = 500)
        numleeches.pred[[i]] <- seq(min(model.data[[i]]$model.data$numleeches, na.rm = TRUE), max(model.data[[i]]$model.data$numleeches, na.rm = TRUE), length.out = 500)

        # empty array to be filled with predictions (some layers will remain unfilled but layers are numbered to correspond with the numbers of the coefficients)
        predC[[i]] <- array(NA, dim = c(500, nsamp[[i]], 6))
        
        # loop over all posterior samples
        for(s in 1:nsamp[[i]]){
            predC[[i]][,s,1] <- plogis(sims[[i]]$mu.eta[s,1] + sims[[i]]$mu.betalpsi1[s] * elev.pred[[i]])   # psi ~ elev
            predC[[i]][,s,5] <- plogis(sims[[i]]$mu.eta[s,1] + sims[[i]]$mu.betalpsi5[s] * reserve.pred[[i]])    # psi ~ reserve
            predC[[i]][,s,6] <- plogis(sims[[i]]$mu.eta[s,2] + sims[[i]]$mu.betalp1[s] * numleeches.pred[[i]])   # p ~ numleeches
        }
        
        # prediction mean
        pmC[[i]] <- apply(predC[[i]], c(1,3), mean)
        # 95% credible interval
        criI[[i]] <- apply(predC[[i]], c(1,3), function(x) quantile(x, prob = c(0.05, 0.95), na.rm=TRUE))

    }

    rm(i)

# draw plots

    pdf("FigS3a_LSU_env_covariates_model23a.pdf", width=3, height=6)
        par(mfrow=c(3,1), mar=c(5,4,1,2))
        # elev
            plot(elev.pred$LSU, pmC$LSU[,1], col="blue", lwd=3, type='l', lty=1, frame=F, ylim=c(0, 0.6), xlab="Elevation", ylab="Community mean occupancy")
            matlines(elev.pred$LSU, t(criI$LSU[,,1]), col="grey", lty=1)
        # reserve
            plot(reserve.pred$LSU, pmC$LSU[,5], col="blue", lwd=3, type='l', lty=1, frame=F, ylim=c(0, 0.8), xlab="Distance to reserve boundary", ylab="Community mean occupancy", yaxt='n')
            axis(side=2, at = seq(0, 0.8, 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8"))
            matlines(reserve.pred$LSU, t(criI$LSU[,,5]), col="grey", lty=1)
        # numleeches
            plot(numleeches.pred$LSU, pmC$LSU[,6], col="blue", lwd=3, type='l', lty=1, frame=F, ylim=c(0, 0.15), xlab="Number of leeches", ylab="Community mean detection", yaxt='n')
            axis(side=2, at = seq(0, 0.15, 0.05), labels = c("0", "0.05", "0.10", "0.15"))
            matlines(numleeches.pred$LSU, t(criI$LSU[,,6]), col="grey", lty=1)
    dev.off()

    pdf("FigS3b_SSU_env_covariates_model23a.pdf", width=3, height=6)
        par(mfrow=c(3,1), mar=c(5,4,1,2))
        # elev
            plot(elev.pred$SSU, pmC$SSU[,1], col="blue", lwd=3, type='l', lty=1, frame=F, ylim=c(0, 0.6), xlab="Elevation", ylab="Community mean occupancy")
            matlines(elev.pred$SSU, t(criI$SSU[,,1]), col="grey", lty=1)
        # reserve
            plot(reserve.pred$SSU, pmC$SSU[,5], col="blue", lwd=3, type='l', lty=1, frame=F, ylim=c(0, 0.8), xlab="Distance to reserve boundary", ylab="Community mean occupancy", yaxt='n')
            axis(side=2, at = seq(0, 0.8, 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8"))
            matlines(reserve.pred$SSU, t(criI$SSU[,,5]), col="grey", lty=1)
        # numleeches
            plot(numleeches.pred$SSU, pmC$SSU[,6], col="blue", lwd=3, type='l', lty=1, frame=F, ylim=c(0, 0.15), xlab="Number of leeches", ylab="Community mean detection", yaxt='n')
            axis(side=2, at = seq(0, 0.15, 0.05), labels = c("0", "0.05", "0.10", "0.15"))
            matlines(numleeches.pred$SSU, t(criI$SSU[,,6]), col="grey", lty=1)
    dev.off()

########################################################################################
## check model convergence for best models ##

    # Rhat values should be close to 1 for convergence, and they are

    # Rhat due to Gelman and Rubin 1992 (aka Brooks-Gelman-Rubin statistic; see also Gelman et al 2014)

    # LSU model23a
        LSU_model23a_Rhat <- bind_rows(
            psi = psi.output$LSU$model23a %>% select(Rhat),
            lpsi0 = lpsi0.output$LSU$model23a %>% select(Rhat),
            lp0 = lp0.output$LSU$model23a %>% select(Rhat),
            occ = occ.output$LSU$model23a %>% select(Rhat),
            Nsite = Nsite.output$LSU$model23a %>% select(Rhat),
            z = z.output$LSU$model23a %>% select(Rhat),
            .id = "model.parameter"
        )
        LSU_model23a_Rhat.plot <- LSU_model23a_Rhat %>%
            ggplot() + geom_histogram(aes(x=Rhat), bins = 12) + facet_wrap(~model.parameter, scales = "free")
        ggexport(LSU_model23a_Rhat.plot, filename = "Ailaoshan_occupancy_MCMC_Rhat_LSU.pdf")

    # SSU model23a
        SSU_model23a_Rhat <- bind_rows(
            psi = psi.output$SSU$model23a %>% select(Rhat),
            lpsi0 = lpsi0.output$SSU$model23a %>% select(Rhat),
            lp0 = lp0.output$SSU$model23a %>% select(Rhat),
            occ = occ.output$SSU$model23a %>% select(Rhat),
            Nsite = Nsite.output$SSU$model23a %>% select(Rhat),
            z = z.output$SSU$model23a %>% select(Rhat),
            .id = "model.parameter"
        )
        SSU_model23a_Rhat.plot <- SSU_model23a_Rhat %>%
            ggplot() + geom_histogram(aes(x=Rhat), bins = 12) + facet_wrap(~model.parameter, scales = "free")
        ggexport(SSU_model23a_Rhat.plot, filename = "Ailaoshan_occupancy_MCMC_Rhat_SSU.pdf")
                
    rm(LSU_model23a_Rhat, SSU_model23a_Rhat)

########################################################################################
