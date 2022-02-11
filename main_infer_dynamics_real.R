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
runname <- "variant_compare_plosbiol"

## Assume individuals can stay detectable for up to 35 days
lastday <- 35

## MCMC settings
nchains <- 3
n_temperatures <- 10
mcmcPars_ct_pt <- list("iterations"=50000,"popt"=0.234,"opt_freq"=1000,
                       "thin"=10,"adaptive_period"=30000,"save_block"=1000,
                       "temperature" = seq(1,101,length.out=n_temperatures),
                       "parallel_tempering_iter" = 5,"max_adaptive_period" = 30000, 
                       "adaptiveLeeway" = 0.2, "max_total_iterations" = 50000)

use_pos <- TRUE ## If TRUE, only uses detectable Ct values. If FALSE, uses all Ct values.

###############################################
## 1) DATA SETUP
###############################################
dat <- read_csv("data/pbio.3001216.s001.csv")

ct_dat_sgtf <- dat %>% filter(is.na(`S gene Ct`)) %>% dplyr::select(`ORF1ab Ct`) %>% 
  rename(ct=`ORF1ab Ct`) %>% mutate(t=35) %>% drop_na() %>% filter(ct < 35) %>% mutate(variant="SGTF")

ct_dat <- dat %>% filter(!is.na(`S gene Ct`)) %>% dplyr::select(`ORF1ab Ct`) %>% 
  rename(ct=`ORF1ab Ct`) %>% mutate(t=35)%>% drop_na() %>% filter(ct < 35)%>% mutate(variant="S_pos")

ct_dat_negs <- tibble(t=35,ct=rep(35, (19176 - 616)),variant="S_pos")

ct_dat <- bind_rows(ct_dat,ct_dat_sgtf)


###############################################
## 2) VIROSOLVER SETUP
###############################################
virosolver_pars <- read.csv("pars/partab_exp_model_compare.csv")

virosolver_pars[virosolver_pars$names %in% c("intercept","true_0"),"values"] <- 35
virosolver_pars[virosolver_pars$names %in% c("intercept","true_0"),"lower_start"] <- 33
virosolver_pars[virosolver_pars$names %in% c("intercept","true_0"),"lower_bound"] <- 33

virosolver_pars[virosolver_pars$names == "level_switch","values"] <- 33
virosolver_pars[virosolver_pars$names == "viral_peak","values"] <- 15
virosolver_pars[virosolver_pars$names == "t_switch","values"] <- 10
virosolver_pars[virosolver_pars$names == "prob_detect","values"] <- 0.25
virosolver_pars[virosolver_pars$names == "overall_prob_scale",c("lower_start","lower_bound")] <- 0
virosolver_pars[virosolver_pars$names == "overall_prob_scale",c("upper_start","upper_bound")] <- 1

virosolver_pars[virosolver_pars$names == "tshift","fixed"] <- 1




if(!use_pos){
  virosolver_pars[virosolver_pars$names == "overall_prob","fixed"] <- 0
  virosolver_pars[virosolver_pars$names == "overall_prob_scale","fixed"] <- 0
  virosolver_pars[virosolver_pars$names == "overall_prob","values"] <- 0.03
  virosolver_pars[virosolver_pars$names == "overall_prob_scale","values"] <- 0.03
} else {
  
  virosolver_pars[virosolver_pars$names == "overall_prob","values"] <- 1
  virosolver_pars[virosolver_pars$names == "overall_prob_scale","values"] <- 1
}

pars <- virosolver_pars$values
names(pars) <- virosolver_pars$names

## Pull out the current values for each parameter, and set these as the prior means
means <- virosolver_pars$values
names(means) <- virosolver_pars$names

## Set standard deviations of prior distribution
sds_gp <- c("obs_sd"=1,"viral_peak"=2,
            "wane_rate2"=2,"t_switch"=2, "level_switch"=2,
            "prob_detect"=0.1,
            "incubation"=0.25, "infectious"=0.5)
## Define prior function
prior_func <- function(pars, ...){
  par_names <- names(pars)
  
  ## Viral kinetics parameters
  obs_sd_prior <- dnorm(pars["obs_sd"], means[which(names(means) == "obs_sd")], sds_gp["obs_sd"],log=TRUE)
  viral_peak_prior <- dnorm(pars["viral_peak"], means[which(names(means) == "viral_peak")], sds_gp["viral_peak"],log=TRUE)
  wane_2_prior <- dnorm(pars["wane_rate2"],means[which(names(means) == "wane_rate2")],sds_gp["wane_rate2"],log=TRUE)
  tswitch_prior <- dnorm(pars["t_switch"],means[which(names(means) == "t_switch")],sds_gp["t_switch"],log=TRUE)
  level_prior <- dnorm(pars["level_switch"],means[which(names(means) == "level_switch")],sds_gp["level_switch"],log=TRUE)
  beta1_mean <- means[which(names(means) == "prob_detect")]
  beta1_sd <- sds_gp["prob_detect"]
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  beta_prior <- dbeta(pars["prob_detect"],beta_alpha,beta_beta,log=TRUE)
  
  #########
  gr_prior <- dnorm(pars["beta"],0,sd=0.1,log=TRUE)
  gr_prior2 <- dnorm(pars["beta_alt"],0,sd=0.1,log=TRUE)
  #########
  
  ## Scale priors, log normal centered on 1
  peak_ct_scale_prior <- dlnorm(pars["viral_peak_scale"],log(1),sd=0.5)
  tswitch_scale_prior <- dlnorm(pars["t_switch_scale"],log(1),sd=0.5)
  
  obs_sd_prior + viral_peak_prior + wane_2_prior + tswitch_prior +
    level_prior + beta_prior + gr_prior + gr_prior2 + 
    peak_ct_scale_prior + tswitch_scale_prior
}

## Test that posterior functions work correctly
f1 <- create_posterior_func_compare(parTab=virosolver_pars,data=ct_dat,PRIOR_FUNC = prior_func,
                                   INCIDENCE_FUNC=virosolver::exponential_growth_model,use_pos = use_pos)
f1(virosolver_pars$values)

chains <- NULL
res <- foreach(j=1:nchains,.packages = c("extraDistr","tidyverse","patchwork","virosolver")) %dopar% {
  devtools::load_all("~/Documents/GitHub/lazymcmc")
  if(!file.exists(paste0("chains/",runname))){
    dir.create(paste0("chains/",runname),recursive = TRUE)
  }
  
  tmp_pars <- virosolver_pars
  startTab <- rep(list(tmp_pars),n_temperatures)
  
  for(k in 1:length(startTab)){
    startTab[[k]] <- generate_viable_start_pars(parTab=tmp_pars,
                                                obs_dat=ct_dat,
                                                CREATE_POSTERIOR_FUNC=create_posterior_func_compare,
                                                INCIDENCE_FUNC=virosolver::exponential_growth_model,
                                                PRIOR_FUNC=prior_func,
                                                use_pos=use_pos,
                                                t_dist=NULL)
  }
  covMat <- diag(nrow(startTab[[1]]))
  mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[[1]][startTab[[1]]$fixed==0,])),w=0.8)
  mvrPars <- rep(list(mvrPars), n_temperatures)
  
  output <- run_MCMC(parTab=startTab,
                     data=ct_dat,
                     INCIDENCE_FUNC=virosolver::exponential_growth_model,
                     PRIOR_FUNC=prior_func,
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

p_incidence1 <- plot_virosolver_comparisons(chains, NULL, NULL,NULL,NULL) + 
  coord_cartesian(ylim=c(0,0.01))

p_compare_kinetics1 <- p_compare_estimated_curves(chains,ages=0:50,N=1,nsamp=100,
                                                  true_pars1=NULL,true_pars2=NULL)

model_func <- create_posterior_func(virosolver_pars,ct_dat,NULL,
                                    INCIDENCE_FUNC=virosolver::exponential_growth_model,solve_ver = "model",
                                    use_pos = use_pos)
model_func(pars)

chains_variant <- chains
chains_variant$beta <- chains_variant$beta
chains_variant$viral_peak <- chains_variant$viral_peak * chains_variant$viral_peak_scale
chains_variant$t_switch <- chains_variant$t_switch * chains_variant$t_switch_scale
chains_variant$overall_prob <- chains_variant$overall_prob_scale

plot_distribution_fits(chains,ct_dat %>% filter(variant=="S_pos"),model_func,pos_only=use_pos)
plot_distribution_fits(chains_variant,ct_dat %>% filter(variant=="SGTF"),model_func,pos_only=use_pos)

ggplot(ct_dat) + geom_dotplot(aes(x=variant,y=ct,group=variant),
                              binaxis="y",stackdir="center",binwidth=1,
                              dotsize=0.25) + 
  scale_y_continuous(trans="reverse")
