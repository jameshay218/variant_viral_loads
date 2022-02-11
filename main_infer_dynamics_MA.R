###############################################
## DISTINGUISHING GROWTH RATES AND VIRAL KINETICS FROM VARIANTS
## - Reads in the Ct value data from Ferguson et al. PLOS Biology from Pillar 2 Birmingham testing
## - Uses virosolver to try to estimate the variant-specific growth rate and kinetics parameters
## James Hay
## jhay@hsph.harvard.edu
###############################################
library(virosolver)
library(tidyverse)
library(patchwork)
library(extraDistr)
library(ggpubr)
library(doParallel)
library(lazymcmc)
#devtools::load_all("~/Documents/GitHub/lazymcmc/")
#devtools::load_all("~/Documents/GitHub/virosolver/")

## Where to perform the simulations
HOME_WD <- "~"
HOME_WD <- "~/Documents/GitHub/"
setwd(paste0(HOME_WD,"/variant_viral_loads"))

## Load functions for line list simulation
source("code/linelist_sim_funcs.R")
source("code/plotting.R")
source("code/seir_funcs.R")
source("code/analysis_funcs.R")
source("code/invasion_rates_KISSLER2020.R")
source("code/simulate_symptomatic_population.R")

## Set up cluster for parallel chains
n_clusters <- 9
cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
registerDoParallel(cl)

###############################################
## KEY ASSUPTIONS
###############################################
runname <- "variant_compare_mass"

## Assume individuals can stay detectable for up to 35 days
lastday <- 35

## MCMC settings
nchains <- 3
n_temperatures <- 10
mcmcPars_ct_pt <- c("iterations"=100000,"popt"=0.44,"opt_freq"=2000,
                    "thin"=50,"adaptive_period"=50000,"save_block"=100)

use_pos <- TRUE ## If TRUE, only uses detectable Ct values. If FALSE, uses all Ct values.

###############################################
## 1) DATA SETUP
###############################################
dat <- read_csv("~/Documents/local_data/BWH/BWH_Panther_Cts_20200329-20211031.csv")
ct_dat <- dat %>% filter((coll_week > as.Date("2020-07-01") & coll_week < as.Date('2020-11-01')) |
                           (coll_week > as.Date("2021-07-01") & coll_week < as.Date('2021-11-01'))) %>%
  filter(Ct_min < 40) %>%
  mutate(variant=ifelse(coll_week < as.Date("2021-01-01"),"A","B")) %>%
  rename(ct=Ct_min) %>%
  group_by(variant) %>%
  mutate(t = as.numeric(coll_week - min(coll_week) + 35)) %>% filter(ct < 38)
times <- 0:max(ct_dat$t)

###############################################
## 2) VIROSOLVER SETUP
###############################################
virosolver_pars <- read.csv("pars/partab_gp_model_compare.csv")
virosolver_pars[virosolver_pars$names %in% c("nu","rho"), "values"] <- c(1.5,0.03)
virosolver_pars[virosolver_pars$names %in% c("nu","rho"), "fixed"] <- 1

virosolver_pars[virosolver_pars$names == "sd_mod","fixed"] <- 1
virosolver_pars[virosolver_pars$names == "sd_mod","values"] <- 1
virosolver_pars[virosolver_pars$names == "sd_mod_wane","values"] <- 14

#virosolver_pars[virosolver_pars$names %in% c("intercept","true_0"),"values"] <- 35
#virosolver_pars[virosolver_pars$names %in% c("intercept","true_0"),"lower_start"] <- 33
#virosolver_pars[virosolver_pars$names %in% c("intercept","true_0"),"lower_bound"] <- 33

virosolver_pars[virosolver_pars$names == "level_switch","values"] <- 38
virosolver_pars[virosolver_pars$names == "level_switch","fixed"] <- 1
#virosolver_pars[virosolver_pars$names == "viral_peak","values"] <- 15
#virosolver_pars[virosolver_pars$names == "t_switch","values"] <- 10
#virosolver_pars[virosolver_pars$names == "prob_detect","values"] <- 0.25
virosolver_pars[virosolver_pars$names == "overall_prob_scale",c("lower_start","lower_bound")] <- 0
virosolver_pars[virosolver_pars$names == "overall_prob_scale",c("upper_start","upper_bound")] <- 1

virosolver_pars[virosolver_pars$names == "tshift","fixed"] <- 1

virosolver_pars[virosolver_pars$names == "desired_mode","fixed"] <- 1
virosolver_pars[virosolver_pars$names == "desired_mode","upper_bound"] <- 25
virosolver_pars[virosolver_pars$names == "desired_mode","lower_bound"] <- 0

virosolver_pars[virosolver_pars$names == "t_switch","values"] <- 13

virosolver_pars[virosolver_pars$names == "overall_prob","values"] <- 1
virosolver_pars[virosolver_pars$names == "overall_prob_scale","values"] <- 1

if(!use_pos){
  virosolver_pars[virosolver_pars$names == "overall_prob","fixed"] <- 0
  virosolver_pars[virosolver_pars$names == "overall_prob_scale","fixed"] <- 0
}

virosolver_pars <- bind_rows(virosolver_pars[!(virosolver_pars$names %in% c("prob","prob_alt")),], 
                             do.call("rbind",replicate(length(times), virosolver_pars[virosolver_pars$names == "prob",],simplify=FALSE)),
                             do.call("rbind",replicate(length(times), virosolver_pars[virosolver_pars$names == "prob_alt",],simplify=FALSE)))
pars <- virosolver_pars$values
names(pars) <- virosolver_pars$names

## Means for priors
means <- virosolver_pars$values
names(means) <- virosolver_pars$names

## Set standard deviations of prior distribution
sds_exp <- sds_seir <- sds_gp <- c("beta"=0.25,"R0"=0.6,
                                   "obs_sd"=0.25,"viral_peak"=1,
                                   "wane_rate2"=1,"t_switch"=1.5,"level_switch"=1,
                                   "prob_detect"=0.015,
                                   "incubation"=0.25, "infectious"=0.5,
                                   "rho"=2,"nu"=0.5,"desired_mode"=0.5)
## Define prior function
prior_func_gp <- function(pars, ...){
  par_names <- names(pars)
  
  ## Viral kinetics parameters
  obs_sd_prior <- dnorm(pars["obs_sd"], means[which(names(means) == "obs_sd")], sds_gp["obs_sd"],log=TRUE)
  viral_peak_prior <- dnorm(pars["viral_peak"], means[which(names(means) == "viral_peak")], sds_gp["viral_peak"],log=TRUE)
  wane_2_prior <- dnorm(pars["wane_rate2"],means[which(names(means) == "wane_rate2")],sds_gp["wane_rate2"],log=TRUE)
  tswitch_prior <- dnorm(pars["t_switch"],means[which(names(means) == "t_switch")],sds_gp["t_switch"],log=TRUE)
  level_prior <- dnorm(pars["level_switch"],means[which(names(means) == "level_switch")],sds_gp["level_switch"],log=TRUE)
  mode_prior <- dnorm(pars["desired_mode"],means[which(names(means) == "desired_mode")],sds_gp["desired_mode"],log=TRUE)
  beta1_mean <- means[which(names(means) == "prob_detect")]
  beta1_sd <- sds_gp["prob_detect"]
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  beta_prior <- dbeta(pars["prob_detect"],beta_alpha,beta_beta,log=TRUE)
  
  #########
  ## Gaussian process prior, un-centered version
  k <- pars[which(par_names %in% c("prob","prob_alt"))]
  ## Leave this - correct for uncentered version as per Chapter 14 Statistical Rethinking
  prob_priors <- sum(dnorm(k, 0, 1, log=TRUE))
  nu_prior <- dexp(pars["nu"], 1/means[which(names(means) == "nu")],log=TRUE)
  rho_prior <- dexp(pars["rho"], 1/means[which(names(means) == "rho")],log=TRUE)
  
  #########
  
  ## Scale priors, log normal centered on 1
  peak_ct_scale_prior <- dlnorm(pars["viral_peak_scale"],log(1),sd=0.5)
  tswitch_scale_prior <- dlnorm(pars["t_switch_scale"],log(1),sd=0.25)
  
  obs_sd_prior + viral_peak_prior + wane_2_prior + tswitch_prior +
    level_prior + beta_prior + #mode_prior +
    prob_priors + nu_prior + rho_prior +
    peak_ct_scale_prior + tswitch_scale_prior
}

## Test that posterior functions work correctly
f1 <- create_posterior_func_compare(parTab=virosolver_pars,data=ct_dat,PRIOR_FUNC = prior_func_gp,
                                   INCIDENCE_FUNC=virosolver::gaussian_process_model,use_pos = use_pos)
f1(virosolver_pars$values)

chains <- NULL
res <- foreach(j=1:nchains,.packages = c("extraDistr","tidyverse","patchwork","virosolver","lazymcmc")) %dopar% {
  if(!file.exists(paste0("chains/",runname))){
    dir.create(paste0("chains/",runname),recursive = TRUE)
  }
  
  
  startTab <- generate_viable_start_pars(virosolver_pars,ct_dat,
                                         create_posterior_func_compare,
                                         virosolver::gaussian_process_model,
                                         prior_func_gp,
                                         t_dist=t_dist,
                                         use_pos=use_pos)
  
  
  output <- run_MCMC(parTab=startTab,
                     data=ct_dat,
                     INCIDENCE_FUNC=virosolver::gaussian_process_model,
                     PRIOR_FUNC=prior_func_gp,
                     mcmcPars=mcmcPars_ct_pt,
                     filename=paste0("chains/",runname,"_",j),
                     CREATE_POSTERIOR_FUNC=create_posterior_func_compare,
                     mvrPars=NULL,
                     OPT_TUNING=0.2,
                     use_pos=use_pos,
                     solve_likelihood=TRUE)
  
  ## Read in chain and remove burn in period
  chain <- read.csv(output$file)
  chain <- chain[chain$sampno > mcmcPars_ct_pt["adaptive_period"],]
  chain$sampno <-chain$sampno + max(chain$sampno)*(j-1)
  chains[[j]] <- chain
}
chains <- do.call("bind_rows",res)
colnames(chains) <- c("sampno",virosolver_pars$names,"lnlike")

p_incidence1 <- plot_virosolver_comparisons(chains, NULL, NULL,NULL,NULL,inc_func=virosolver::gaussian_process_model,times=times) #+ 
  #coord_cartesian(ylim=c(0,0.01))

p_compare_kinetics1 <- p_compare_estimated_curves(chains,ages=0:50,N=1,nsamp=100,
                                                  true_pars1=NULL,true_pars2=NULL)

model_func <- create_posterior_func(virosolver_pars,ct_dat,NULL,
                                    INCIDENCE_FUNC=virosolver::gaussian_process_model,solve_ver = "model",
                                    use_pos = use_pos)
model_func(pars)

chains_variant <- chains
chains_variant$beta <- chains_variant$beta
chains_variant[,which(colnames(chains_variant)=="prob_alt")] <- chains_variant[,which(colnames(chains_variant)=="prob")] 
chains_variant$viral_peak <- chains_variant$viral_peak * chains_variant$viral_peak_scale
chains_variant$t_switch <- chains_variant$t_switch * chains_variant$t_switch_scale
chains_variant$overall_prob <- chains_variant$overall_prob_scale

plot_distribution_fits(chains,ct_dat %>% filter(variant=="A"),model_func,pos_only=use_pos)
plot_distribution_fits(chains_variant,ct_dat %>% filter(variant=="B"),model_func,pos_only=use_pos)

ggplot(ct_dat) + geom_dotplot(aes(x=variant,y=ct,group=variant),
                              binaxis="y",stackdir="center",binwidth=1,
                              dotsize=0.25) + 
  scale_y_continuous(trans="reverse")
