###############################################
## DISTINGUISHING GROWTH RATES AND VIRAL KINETICS FROM VARIANTS
## - Simulates Ct values for two variants with different viral kinetics and growth rates
## - Simulate samples taken under random cross-sectional surveillance early on when two viruses
## - are seeded at the same time, or later on if one variant is seeded later than the other.
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
## Assume individuals can stay detectable for up to 35 days
lastday <- 35

## Day of the simulation where the early and late samples are taken
samp_time_early <- 50
samp_time_late <- 270

## How much lower is the peak Ct of the new variant?
ct_diff <- 5
## How many more days does it take for the new variant to be cleared?
tswitch_diff <- 5

## How many Ct values to simulate at each sample time?
sample_size <- 100

## MCMC settings
nchains <- 3
n_temperatures <- 10
mcmcPars_ct_pt <- list("iterations"=50000,"popt"=0.234,"opt_freq"=1000,
                       "thin"=10,"adaptive_period"=50000,"save_block"=1000,
                       "temperature" = seq(1,101,length.out=n_temperatures),
                       "parallel_tempering_iter" = 5,"max_adaptive_period" = 50000, 
                       "adaptiveLeeway" = 0.2, "max_total_iterations" = 50000)

use_pos <- TRUE ## If TRUE, only uses detectable Ct values. If FALSE, uses all Ct values.

###############################################
## 1) SEIR SIMULATION
###############################################
## SEIR parameters
pars_same_seed <- c(sigma1.val = 0,#1/(45*7*2), ## Immune waning to strain 1
          sigma2.val = 0,#1/(45*7*2), ## Immune waning to strain 2
          nu.val = 1/3, ## Latent period
          gamma.val = 1/7, ## Infectious period
          chi12.val = 0.75, ## Cross immunity strain 1 confers against strain 2
          chi21.val = 0.75, ## Cross immunity strain 2 confers against strain 1
          beta.val1=1.5/7, ## Transmission rate of strain 1 (R0*gamma)
          beta.val2=2.5/7, ## Transmission rate of strain 2 (R0*gamma)
          kappa.val = 1/100000, ## Daily importation rate of infected individuals
          importtime1 = 0, ## Time of importation of strain 1
          importtime2 = 0, ## Time of importation of strain 2
          importlength = 7) ## Duration of importations

## Parameters with new variant being seeded at the same time or later on
pars_late_seed <- pars_same_seed
pars_late_seed["importtime2"] <- 180

states <- c(S1S2 = 1,E1S2 = 0,S1E2 = 0,E1E2 = 0,I1S2 = 0, 
            S1I2 = 0, R1S2 = 0,I1E2 = 0, E1I2 = 0, S1R2 = 0, 
            R1E2 = 0, I1I2 = 0, E1R2 = 0, R1I2 = 0, I1R2 = 0, 
            R1R2 = 0, inc1 = 0, inc2 = 0) # Initial conditions

times <- seq(0, 365*1.5,by=1) ## Run model for 1.5 years


########################################################
## A) SEIR SIMULATION if seeded at the same time 
## Solve SEIR model and pull out variant-specific incidence curves
seir_dynamics_same <- run_2strain_seir_simulation(pars_same_seed, states,times)
virus1_inc_same <- seir_dynamics_same$virus1_inc
virus2_inc_same <- seir_dynamics_same$virus2_inc
virus_inc_same <- virus1_inc_same + virus2_inc_same

## Calculate daily growth rates of the viruses
gr1_same <- tibble(t=times,gr=c(0,log(virus1_inc_same[2:length(virus1_inc_same)]/virus1_inc_same[1:(length(virus1_inc_same)-1)])),inc=virus1_inc_same,virus="Original variant")
gr2_same <- tibble(t=times,gr=c(0,log(virus2_inc_same[2:length(virus2_inc_same)]/virus2_inc_same[1:(length(virus2_inc_same)-1)])),inc=virus2_inc_same,virus="New variant")
gr_overall_same <- tibble(t=times,gr=c(0,log(virus_inc_same[2:length(virus_inc_same)]/virus_inc_same[1:(length(virus_inc_same)-1)])),inc=virus_inc_same,virus="Overall")
grs_same <- bind_rows(gr1_same, gr2_same,gr_overall_same)

## Find the average 35-day growth rate in the scenario with seeding on the same day
tmp1_same <- (virus1_inc_same[seq(samp_time_early-lastday,samp_time_early,by=1)])
mean_gr_virus1_same <- mean(log(tmp1_same[2:length(tmp1_same)]/tmp1_same[1:(length(tmp1_same)-1)]))
tmp2_same <- (virus2_inc_same[seq(samp_time_early-lastday,samp_time_early,by=1)])
mean_gr_virus2_same <- mean(log(tmp2_same[2:length(tmp2_same)]/tmp2_same[1:(length(tmp2_same)-1)]))

## Put true values into tibble for later plotting
true_vals_same <- tibble(gr=c(mean_gr_virus1_same,mean_gr_virus2_same,mean_gr_virus2_same-mean_gr_virus1_same),
                         name=c("Original variant","New variant","Difference")) 
true_vals_same$name <- factor(true_vals_same$name,levels=c("Original variant","New variant","Difference"))
real_v1_gr_same <- tibble(t=0:lastday,prob_infection=tmp1_same)
real_v2_gr_same <- tibble(t=0:lastday,prob_infection=tmp2_same)

if(use_pos){
  real_v1_gr_same$prob_infection <- real_v1_gr_same$prob_infection/sum(real_v1_gr_same$prob_infection)
  real_v2_gr_same$prob_infection <- real_v2_gr_same$prob_infection/sum(real_v2_gr_same$prob_infection)
}

########################################################
## B) SEIR SIMULATION if seeded at a later time
## Solve SEIR model and pull out variant-specific incidence curves
seir_dynamics_late <- run_2strain_seir_simulation(pars_late_seed, states,times)
virus1_inc_late <- seir_dynamics_late$virus1_inc
virus2_inc_late <- seir_dynamics_late$virus2_inc
virus_inc_late <- virus1_inc_late + virus2_inc_late

## Calculate daily growth rates of the viruses
gr1_late <- tibble(t=times,gr=c(0,log(virus1_inc_late[2:length(virus1_inc_late)]/virus1_inc_late[1:(length(virus1_inc_late)-1)])),inc=virus1_inc_late,virus="Original variant")
gr2_late <- tibble(t=times,gr=c(0,log(virus2_inc_late[2:length(virus2_inc_late)]/virus2_inc_late[1:(length(virus2_inc_late)-1)])),inc=virus2_inc_late,virus="New variant")
gr_overall_late <- tibble(t=times,gr=c(0,log(virus_inc_late[2:length(virus_inc_late)]/virus_inc_late[1:(length(virus_inc_late)-1)])),inc=virus_inc_late,virus="Overall")
grs_late <- bind_rows(gr1_late, gr2_late,gr_overall_late)

## Find the average 35-day growth rate in the scenario with seeding on the same day
tmp1_late <- (virus1_inc_late[seq(samp_time_late-lastday,samp_time_late,by=1)])
mean_gr_virus1_late <- mean(log(tmp1_late[2:length(tmp1_late)]/tmp1_late[1:(length(tmp1_late)-1)]))
tmp2_late <- (virus2_inc_late[seq(samp_time_late-lastday,samp_time_late,by=1)])
mean_gr_virus2_late <- mean(log(tmp2_late[2:length(tmp2_late)]/tmp2_late[1:(length(tmp2_late)-1)]))

## Put true values into tibble for later plotting
true_vals_late <- tibble(gr=c(mean_gr_virus1_late,mean_gr_virus2_late,mean_gr_virus2_late-mean_gr_virus1_late),
                         name=c("Original variant","New variant","Difference")) 
true_vals_late$name <- factor(true_vals_late$name,levels=c("Original variant","New variant","Difference"))
real_v1_gr_late <- tibble(t=0:lastday,prob_infection=tmp1_late/sum(tmp1_late))
real_v2_gr_late <- tibble(t=0:lastday,prob_infection=tmp2_late/sum(tmp2_late))


###############################################
## 2) VIRAL KINETICS
###############################################
ages <- 1:lastday

model_pars <- read.csv("pars/partab_seir_model.csv")
vl_pars <- model_pars$values
names(vl_pars) <- model_pars$names

## Have version with both persistent and lower Cts
vl_pars_both <- vl_pars
vl_pars_both["t_switch"] <- vl_pars_both["t_switch"] + tswitch_diff
vl_pars_both["viral_peak"] <- vl_pars_both["viral_peak"] - ct_diff

vl_pars1 <- vl_pars
vl_pars2 <- vl_pars_both

## Tibble for later plotting with true ratio of viral kinetics parameters
scale_key <- c("viral_peak_scale"="Relative peak\n Ct value","t_switch_scale"="Relative duration\n of initial waning")
real_scales <- tibble(name=scale_key,value=c(vl_pars2["viral_peak"]/vl_pars1["viral_peak"], vl_pars2["t_switch"]/vl_pars1["t_switch"]))

########################################################
## A) SIMULATE CT VALUES if seeded at the same time
cts_1_same <- tibble(ct=simulate_cross_section(vl_pars1, ages, virus1_inc_same,obs_time=samp_time_early,N=sample_size,use_pos=use_pos),variant="Original variant")
cts_2_same <- tibble(ct=simulate_cross_section(vl_pars2, ages, virus2_inc_same,obs_time=samp_time_early,N=sample_size,use_pos=use_pos),variant="New variant")
cts_sim_comb_same <- bind_rows(cts_1_same,cts_2_same) %>% mutate(t=lastday)
cts_sim_comb_same$virus <- factor(cts_sim_comb_same$variant,levels=variant_levels)

## Dotplot of the simulated Ct values
## Define a temporary function
TMP_dotplot_function <- function(ct_data){
  p <- ggplot(ct_data) + 
    geom_violin(aes(x=virus,y=ct,fill=virus),alpha=0.1,trim=FALSE,col=NA,width=0.75) +
    geom_dotplot(aes(x=virus,y=ct,fill=virus),binaxis="y",stackdir="center",binwidth=1,dotsize=1) + 
    geom_errorbar(data=cts_sim_comb_same%>%group_by(virus)%>%summarize(median_ct=median(ct)),
                  aes(y=median_ct,ymin=median_ct,ymax=median_ct,x=virus),size=0.5,col="grey10") +
    variant_fill_scale + variant_color_scale +
    scale_y_continuous(trans="reverse",limits=c(40,7)) +
    ylab("Ct value") +
    theme_overall + theme_nice_axes + theme(legend.position="none",axis.title.x=element_blank(),plot.title=element_text(size=7)) +
    scale_x_discrete(labels = c("Original variant" = "Original variant",
                                "New variant"="New variant")) 
}
## Call this function
p_ct_samp_same <- TMP_dotplot_function(cts_sim_comb_same)

## B) SIMULATE CT VALUES if seeded at a later time
cts_1_late <- tibble(ct=simulate_cross_section(vl_pars1, ages, virus1_inc_late,obs_time=samp_time_late,N=sample_size,use_pos=use_pos),variant="Original variant")
cts_2_late <- tibble(ct=simulate_cross_section(vl_pars2, ages, virus2_inc_late,obs_time=samp_time_late,N=sample_size,use_pos=use_pos),variant="New variant")
cts_sim_comb_late <- bind_rows(cts_1_late,cts_2_late) %>% mutate(t=35)
cts_sim_comb_late$virus <- factor(cts_sim_comb_late$variant,levels=variant_levels)

## Dotplot of the simulated Ct values
p_ct_samp_late <- TMP_dotplot_function(cts_sim_comb_late)

###############################################
## 3) VIROSOLVER SETUP
###############################################
virosolver_pars <- read.csv("pars/partab_exp_model_compare.csv")

if(!use_pos){
  virosolver_pars[virosolver_pars$names == "overall_prob","fixed"] <- 0
}

pars <- virosolver_pars$values
names(pars) <- virosolver_pars$names

## Pull out the current values for each parameter, and set these as the prior means
means <- virosolver_pars$values
names(means) <- virosolver_pars$names

## Set standard deviations of prior distribution
sds_gp <- c("obs_sd"=0.5,"viral_peak"=2,
            "wane_rate2"=1,"t_switch"=3,"level_switch"=1,
            "prob_detect"=0.03,
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
f1 <- create_posterior_func_compare(parTab=virosolver_pars,data=cts_sim_comb_same,PRIOR_FUNC = prior_func,
                                   INCIDENCE_FUNC=virosolver::exponential_growth_model,use_pos = use_pos)
f1(virosolver_pars$values)
f2 <- create_posterior_func_compare(parTab=virosolver_pars,data=cts_sim_comb_late,PRIOR_FUNC = prior_func,
                                   INCIDENCE_FUNC=virosolver::exponential_growth_model,use_pos = use_pos)
f2(virosolver_pars$values)


data_list <- rep(list(cts_sim_comb_same,cts_sim_comb_late,cts_sim_comb_same),each=nchains)
runnames <- rep(c("virosolver_same", "virosolver_late","virosolver_same_fixed_pars"),each=nchains)

if(!use_pos){
  runnames <- paste0(runnames, "all_cts")
}

chain_nos <- rep(1:nchains,3)

chains <- NULL
res <- foreach(j=seq_along(chain_nos),.packages = c("extraDistr","tidyverse","patchwork","virosolver","lazymcmc")) %dopar% {
  data_use <- data_list[[j]]
  chain_no <- chain_nos[j]
  runname <- runnames[j]
  
  if(!file.exists(paste0("chains/",runname))){
    dir.create(paste0("chains/",runname),recursive = TRUE)
  }
  
  tmp_pars <- virosolver_pars
  if(runname == "virosolver_same_fixed_pars"){
    tmp_pars[!(tmp_pars$names %in% c("beta","beta_alt","viral_peak_scale","t_switch_scale")),"fixed"] <- 1
  }
  startTab <- rep(list(tmp_pars),n_temperatures)
  for(k in 1:length(startTab)){
    startTab[[k]] <- generate_viable_start_pars(parTab=tmp_pars,
                                                obs_dat=data_use,
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
                     data=data_use,
                     INCIDENCE_FUNC=virosolver::exponential_growth_model,
                     PRIOR_FUNC=prior_func,
                     mcmcPars=mcmcPars_ct_pt,
                     filename=paste0("chains/",runname,"/",runname,"_",chain_no),
                     CREATE_POSTERIOR_FUNC=create_posterior_func_compare,
                     mvrPars=mvrPars,
                     OPT_TUNING=0.2,
                     use_pos=use_pos,
                     solve_likelihood=TRUE)
  
  ## Read in chain and remove burn in period
  chain <- read.csv(output$file)
  chain <- chain[chain$sampno > mcmcPars_ct_pt["adaptive_period"],]
  chain$sampno <-chain$sampno + max(chain$sampno)*(chain_no-1)
  chains[[j]] <- chain
}
chains_same <- do.call("bind_rows",res[1:nchains])
chains_late <- do.call("bind_rows",res[(nchains+1):(nchains+nchains)])
chains_fixed <- do.call("bind_rows",res[(nchains+nchains+1):(3*nchains)])


p_incidence1 <- plot_virosolver_comparisons(chains_same, true_vals_same, real_scales,real_v1_gr_same,real_v2_gr_same) + 
  coord_cartesian(ylim=c(0,0.15))
p_incidence2 <- plot_virosolver_comparisons(chains_late, true_vals_late, real_scales,real_v1_gr_late,real_v2_gr_late)+ 
  coord_cartesian(ylim=c(0,0.15))
p_incidence3 <- plot_virosolver_comparisons(chains_fixed, true_vals_same, real_scales,real_v1_gr_same,real_v2_gr_same)

p_compare_kinetics1 <- p_compare_estimated_curves(chains_same,ages=0:lastday,N=100,nsamp=1000,true_pars1=vl_pars1,true_pars2=vl_pars2)
p_compare_kinetics2 <- p_compare_estimated_curves(chains_late,ages=0:lastday,N=100,nsamp=1000,true_pars1=vl_pars1,true_pars2=vl_pars2)
p_compare_kinetics3 <- p_compare_estimated_curves(chains_late,ages=0:lastday,N=100,nsamp=1000,true_pars1=vl_pars1,true_pars2=vl_pars2)

p_RHS <- (p_ct_samp_same+labs(tag="D")+ggtitle("Both increasing, new variant at higher rate")+theme(plot.title = element_text(size=8,face="bold"))) / 
  (p_incidence1+labs(tag="E")) / (p_compare_kinetics1+labs(tag="F") + theme(legend.position="none"))
p_LHS <- (p_ct_samp_late+labs(tag="A")+ggtitle("Original variant decreasing, new variant increasing")+theme(plot.title = element_text(size=8,face="bold"))) / 
  (p_incidence2+labs(tag="B")) / (p_compare_kinetics2+labs(tag="C")+ theme(legend.position="none"))


p_main <- p_LHS | p_RHS
if(FALSE){
  ggsave("figures/estimate_diffs.png",p_main,height=8,width=8,units="in",dpi=300)
  ggsave("figures/estimate_diffs.pdf",p_main,height=8,width=8)
}

