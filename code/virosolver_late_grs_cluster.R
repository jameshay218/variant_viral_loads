library(virosolver)
library(tidyverse)
library(patchwork)
library(extraDistr)
library(ggpubr)
library(doParallel)

## Where to perform the simulations
HOME_WD <- "~"
#HOME_WD <- "~/Documents/GitHub/"
setwd(HOME_WD)

devtools::load_all(paste0(HOME_WD,"/lazymcmc"))

## Load functions for line list simulation
source(paste0(HOME_WD,"/variant_viral_loads/code/linelist_sim_funcs.R"))
source(paste0(HOME_WD,"/variant_viral_loads/code/plotting.R"))
source(paste0(HOME_WD,"/variant_viral_loads/code/seir_funcs.R"))
source(paste0(HOME_WD,"/variant_viral_loads/code/analysis_funcs.R"))
source(paste0(HOME_WD,"/variant_viral_loads/code/invasion_rates_KISSLER2020.R"))
source(paste0(HOME_WD,"/variant_viral_loads/code/simulate_symptomatic_population.R"))

## Creating and Setting Directories:
main_wd <- paste0(HOME_WD,"/variant_viral_loads/")
chainwd <- paste0(HOME_WD, "/variant_viral_loads/mcmc_chains/virosolver_late")
plot_wd <- paste0(HOME_WD, "/variant_viral_loads/plots/virosolver_late")
data_wd <- paste0(HOME_WD, "/variant_viral_loads/data/virosolver_late")
results_wd <- paste0(HOME_WD, "/variant_viral_loads/results/virosolver_late")

if(!file.exists(chainwd)) dir.create(chainwd,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)
if(!file.exists(data_wd)) dir.create(data_wd,recursive = TRUE)
if(!file.exists(results_wd)) dir.create(results_wd,recursive = TRUE)

## Where to perform the simulations
setwd(main_wd)

n_clusters <- 3
cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
registerDoParallel(cl)


########################################
## 2. Simulation settings
########################################
## Arguments for this run. Can generate automatically here
control_table <- expand.grid(repeat_no=1:100,samp_size=c(25,50,100,250,500))
control_table$run_name <- 1:nrow(control_table)

n_temperatures <- 10
mcmcPars_ct_pt <- list("iterations"=5000,"popt"=0.234,"opt_freq"=1000,
                       "thin"=1,"adaptive_period"=5000,"save_block"=1000,
                       "temperature" = seq(1,101,length.out=n_temperatures),
                       "parallel_tempering_iter" = 5,"max_adaptive_period" = 5000, 
                       "adaptiveLeeway" = 0.2, "max_total_iterations" = 5000)

nchains <- 3
n_samp <- 100 ## Number of posterior samples for plots


## Set Simulation Number and get sim settings
Sim <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#Sim <- 1
print(paste0("Starting Simulation Number: ",Sim))

## Name of this run
runname <- control_table$run_name[Sim] 
## Which simulated population to use?
repeat_no <- control_table$repeat_no[Sim] 
## How many Ct values for each variant to sample?
samp_size <- control_table$samp_size[Sim] 

set.seed(Sim)

runname_use <- paste0("run_",samp_size,"_",runname)
dir.create(paste0(chainwd,"/",runname_use),recursive = TRUE)

########################################
## 2. SEIR simulation
########################################
samp_time <- 270

## SEIR parameters
pars <- c(sigma1.val = 0,#1/(45*7*2), ## Immune waning to strain 1
          sigma2.val = 0,#1/(45*7*2), ## Immune waning to strain 2
          nu.val = 1/3, ## Latent period
          gamma.val = 1/7, ## Infectious period
          chi12.val = 0.75, ## Cross immunity strain 1 confers against strain 2
          chi21.val = 0.75, ## Cross immunity strain 2 confers against strain 1
          #amplitude = 0, 
          #baseline = 3, 
          #phi.val = 0, 
          beta.val1=1.5/7, ## Transmission rate of strain 1 (R0*gamma)
          beta.val2=2.5/7, ## Transmission rate of strain 2 (R0*gamma)
          kappa.val = 1/100000, ## Daily importation rate of infected individuals
          importtime1 = 0, ## Time of importation of strain 1
          importtime2 = 180, ## Time of importation of strain 2
          importlength = 7) ## Duration of importations

states <- c(S1S2 = 1,E1S2 = 0,S1E2 = 0,E1E2 = 0,I1S2 = 0, 
            S1I2 = 0, R1S2 = 0,I1E2 = 0, E1I2 = 0, S1R2 = 0, 
            R1E2 = 0, I1I2 = 0, E1R2 = 0, R1I2 = 0, I1R2 = 0, 
            R1R2 = 0, inc1 = 0, inc2 = 0) # Initial conditions

times <- seq(0, 365*1.5,by=1) ## Run model for 1.5 years

seir_dynamics <- run_2strain_seir_simulation(pars, states,times)
virus1_inc <- seir_dynamics$virus1_inc
virus2_inc <- seir_dynamics$virus2_inc
virus_inc <- virus1_inc + virus2_inc


tmp1 <- (virus1_inc[seq(samp_time-35,samp_time,by=1)])
mean_gr_virus1 <- mean(log(tmp1[2:length(tmp1)]/tmp1[1:(length(tmp1)-1)]))
tmp2 <- (virus2_inc[seq(samp_time-35,samp_time,by=1)])
mean_gr_virus2 <- mean(log(tmp2[2:length(tmp2)]/tmp2[1:(length(tmp2)-1)]))

########################################
## 4. Simulate Ct values
########################################
## Assume individuals can stay detectable for up to 35 days
lastday <- 35
ages <- 1:lastday

model_pars <- read.csv("pars/partab_seir_model.csv")
vl_pars <- model_pars$values
names(vl_pars) <- model_pars$names
vl_pars["tshift"] <- 1
vl_pars["true_0"] <- 45
vl_pars["desired_mode"] <- 4
## Have version with lower peak Ct
vl_pars_peak <- vl_pars
vl_pars_peak["viral_peak"] <- vl_pars_peak["viral_peak"] - 5

## Have version with more persistent Cts
vl_pars_persistent <- vl_pars
vl_pars_persistent["t_switch"] <- vl_pars_persistent["t_switch"] + 5

## Have version with both persistent and lower Cts
vl_pars_both <- vl_pars
vl_pars_both["t_switch"] <- vl_pars_both["t_switch"] + 5
vl_pars_both["viral_peak"] <- vl_pars_both["viral_peak"] - 5

vl_pars1 <- vl_pars
vl_pars2 <- vl_pars_both
cts_1 <- tibble(ct=simulate_cross_section(vl_pars1, ages, virus1_inc,obs_time=samp_time,N=samp_size),variant="Original variant")
cts_2 <- tibble(ct=simulate_cross_section(vl_pars2, ages, virus2_inc,obs_time=samp_time,N=samp_size),variant="New variant, different kinetics")
cts_sim_comb <- bind_rows(cts_1,cts_2) %>% mutate(t=35)
cts_sim_comb$virus <- factor(cts_sim_comb$variant,levels=variant_levels)

########################################
## 5. MCMC settings
########################################
virosolver_pars <- read.csv(paste0(main_wd,"pars/partab_exp_model_compare.csv"))

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
  
  ## Scale priors
  peak_ct_scale_prior <- dlnorm(pars["viral_peak_scale"],log(1),sd=0.5)
  tswitch_scale_prior <- dlnorm(pars["t_switch_scale"],log(1),sd=0.5)
  
  
  obs_sd_prior + viral_peak_prior + wane_2_prior + tswitch_prior +
    level_prior + beta_prior + gr_prior + gr_prior2 + 
    peak_ct_scale_prior + tswitch_scale_prior
}


f <- create_posterior_func_compare(parTab=virosolver_pars,data=cts_sim_comb,PRIOR_FUNC = prior_func,INCIDENCE_FUNC=virosolver::exponential_growth_model,use_pos = TRUE)
f(virosolver_pars$values)


########################################
## 6. Run MCMC
########################################
## Run for each chain
chains <- NULL
res <- foreach(j=1:nchains,.packages = c("extraDistr","tidyverse","patchwork","virosolver")) %dopar% {
  devtools::load_all(paste0(HOME_WD,"/lazymcmc"))
  
  startTab <- rep(list(virosolver_pars),n_temperatures)
  for(k in 1:length(startTab)){
    startTab[[k]] <- generate_viable_start_pars(parTab=virosolver_pars,
                                                obs_dat=cts_sim_comb,
                                                CREATE_POSTERIOR_FUNC=create_posterior_func_compare,
                                                INCIDENCE_FUNC=virosolver::exponential_growth_model,
                                                PRIOR_FUNC=prior_func,
                                                use_pos=TRUE,
                                                t_dist=NULL)
  }
  covMat <- diag(nrow(startTab[[1]]))
  mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[[1]][startTab[[1]]$fixed==0,])),w=0.8)
  mvrPars <- rep(list(mvrPars), n_temperatures)
  
  
  output <- run_MCMC(parTab=startTab,
                     data=cts_sim_comb,
                     INCIDENCE_FUNC=virosolver::exponential_growth_model,
                     PRIOR_FUNC=prior_func,
                     mcmcPars=mcmcPars_ct_pt,
                     filename=paste0(chainwd,"/",runname_use,"/",runname_use,"_",j),
                     CREATE_POSTERIOR_FUNC=create_posterior_func_compare,
                     mvrPars=mvrPars,
                     OPT_TUNING=0.2,
                     use_pos=TRUE,
                     solve_likelihood=TRUE)
  ## Read in chain and remove burn in period
  chain <- read.csv(output$file)
  chain <- chain[chain$sampno > mcmcPars_ct_pt["adaptive_period"],]
  chain$sampno <-chain$sampno + max(chain$sampno)*(j-1)
  chains[[j]] <- chain
}
chain <- do.call("bind_rows",res)

chain$beta_diff <- chain$beta_alt - chain$beta
chain_grs <- chain %>% dplyr::select(sampno, beta, beta_alt,beta_diff) %>% pivot_longer(-sampno)
gr_key <- c("beta"="Original variant","beta_alt"="New variant", "beta_diff"="Difference")
chain_grs$name <- gr_key[chain_grs$name]
chain_grs$name <- factor(chain_grs$name,levels=c("Original variant","New variant","Difference"))
true_vals <- tibble(gr=c(mean_gr_virus1,mean_gr_virus2,mean_gr_virus2-mean_gr_virus1),name=c("Original variant","New variant","Difference")) 
true_vals$name <- factor(true_vals$name,levels=c("Original variant","New variant","Difference"))

p1 <- ggplot(chain_grs) + 
  geom_hline(yintercept=0,linetype="dashed")+
  geom_violin(aes(x=name,y=value,fill=name),alpha=0.5,draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_point(data=true_vals,aes(x=name,y=gr,col=name)) +
  variant_fill_scale + variant_color_scale +
  theme_overall +
  scale_y_continuous(limits=c(-0.25,0.25),breaks=seq(-0.5,0.5,by=0.1)) +
  ylab("Growth rate estimate") +
  xlab("")


chain_grs$run_ID <- runname
chain_grs$samp_size <- samp_size


chain_scales <- chain %>% dplyr::select(sampno, viral_peak_scale, t_switch_scale) %>% pivot_longer(-sampno)
scale_key <- c("viral_peak_scale"="Relative peak\n Ct value","t_switch_scale"="Relative duration\n of initial waning")
chain_scales$name <- scale_key[chain_scales$name]
chain_scales$run_ID <- runname
chain_scales$samp_size <- samp_size


real_v1_gr <- tibble(t=0:35,prob_infection=tmp1/sum(tmp1))
v1_preds <- virosolver::plot_prob_infection(chain, 100,exponential_growth_model, 0:35,true_prob_infection = real_v1_gr)
v1_preds_quants <- v1_preds$predictions %>% group_by(t) %>% summarize(lower95=quantile(prob_infection,0.025),lower50=quantile(prob_infection,0.25),
                                                                      median=median(prob_infection),
                                                                      upper50=quantile(prob_infection,0.75),upper95=quantile(prob_infection,0.975)) 
chain1 <- chain
chain1$beta <- chain1$beta_alt
real_v2_gr <- tibble(t=0:35,prob_infection=tmp2/sum(tmp2))
v2_preds <- virosolver::plot_prob_infection(chain1, 100,exponential_growth_model, 0:35,true_prob_infection = real_v2_gr)
v2_preds_quants <- v2_preds$predictions %>% group_by(t) %>% summarize(lower95=quantile(prob_infection,0.025),lower50=quantile(prob_infection,0.25),
                                                                      median=median(prob_infection),
                                                                      upper50=quantile(prob_infection,0.75),upper95=quantile(prob_infection,0.975)) 
v1_preds_quants$virus <- "Original variant"
v2_preds_quants$virus <- "New variant, different kinetics"

v_preds_comb <- bind_rows(v1_preds_quants,v2_preds_quants)
v_preds_comb$run_ID <- runname
v_preds_comb$samp_size <- samp_size


if(!file.exists(paste0(results_wd,"/incidence_curves/"))) dir.create(paste0(results_wd,"/incidence_curves/"),recursive = TRUE)
if(!file.exists(paste0(results_wd,"/viral_kinetics/"))) dir.create(paste0(results_wd,"/viral_kinetics/"),recursive = TRUE)
if(!file.exists(paste0(results_wd,"/growth_rates/"))) dir.create(paste0(results_wd,"/growth_rates/"),recursive = TRUE)

write.csv(v_preds_comb, paste0(results_wd,"/incidence_curves/incidence_curves_",samp_size,"_",runname,".csv"))
write.csv(chain_scales, paste0(results_wd,"/viral_kinetics/viral_kinetics_",samp_size,"_",runname,".csv"))
write.csv(chain_grs, paste0(results_wd,"/growth_rates/growth_rates_",samp_size,"_",runname,".csv"))
