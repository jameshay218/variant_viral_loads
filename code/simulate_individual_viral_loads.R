#################################################################
## SIMULATE MANY SWAB VIRAL LOAD TRAJECTORIES WITH AND WITHOUT OBSERVATION ERROR
#################################################################
## 1. Reads in the MCMC fits to the Wolfel et al. SWAB data
## 2. Simulates viral load trajectories with observation error for each individual
library(tidyverse)
library(patchwork)
library(extraDistr)
library(virosolver)
library(ggpubr)
library(lazymcmc)
library(rethinking)

HOME_WD <- "~"
HOME_WD <- "~/Documents/GitHub/variant_viral_loads/"
setwd(HOME_WD)
source("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/functions/simulation_functions.R")
source("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/functions/model_funcs_multivariate_hinge.R")

set.seed(0)

## Sample size
samp_size <- 100000
## Duration of epidemic in days
times <- 0:365
run_name <- "swab"
shift_twane <- 0 ## Can make the waning duration average shift_twane days longer
savewd <- "~/Google Drive/nCoV/sims_for_brian/"

## Viral kinetics pars
## Change wd to local path
parTab <- read.csv("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/pars/partab_multivariate_hinge.csv",stringsAsFactors=FALSE)
chains <- load_mcmc_chains(paste0("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/chains/chains_swab"), parTab, 
                           FALSE, 100, burnin=1000000,multi=TRUE)
chain <- as.data.frame(chains[[2]])
chain$sampno <- 1:nrow(chain)
chain$wane_mean <- chain$wane_mean + shift_twane

chain_vl2 <- chain
chain_vl2$viral_peak_mean <- chain_vl2$viral_peak_mean + 1
chain_vl2$wane_mean <- chain_vl2$wane_mean + 5

cts_indiv_vl1 <- simulate_individual_level_data(samp_size, virus1_inc, times, chain, parTab, max_vl=16)
cts_indiv_vl1 <- cts_indiv_vl1 %>% mutate(virus="Original variant")
cts_indiv_vl2 <- simulate_individual_level_data(samp_size, virus2_inc, times, chain, parTab, max_vl=16)
cts_indiv_vl2 <- cts_indiv_vl2 %>% mutate(virus="New variant, same kinetics") %>% mutate(i = i + max(cts_indiv_vl1$i))
cts_indiv_vl2_alt <- simulate_individual_level_data(samp_size, virus2_inc, times, chain_vl2, parTab, max_vl=16)
cts_indiv_vl2_alt <- cts_indiv_vl2_alt %>% mutate(virus="New variant, different kinetics") %>% mutate(i = i + max(cts_indiv_vl2$i))

cts_indiv_comb <- bind_rows(cts_indiv_vl1,cts_indiv_vl2_alt) %>% bind_rows(cts_indiv_vl2) %>% mutate(virus=factor(virus,levels=variant_levels))
cts_indiv_sub <- cts_indiv_comb %>% filter(days_since_infection >= 0, days_since_infection <= 35, obs >= 0) %>% filter(ct < 40)

obs_delays <- tibble(i=unique(cts_indiv_comb$i),confirm_delay=extraDistr::rdgamma(length(unique(cts_indiv_comb$i)),5, 0.8))

hist(obs_delays$confirm_delay)

cts_indiv_sub <- cts_indiv_comb %>% left_join(obs_delays) %>% 
  mutate(obs_time = round(inf_time+incu_period+confirm_delay)) %>% 
  mutate(onset_time=inf_time+incu_period) %>%
  mutate(days_since_onset=obs_time-onset_time) 
cts_indiv_sub %>% 
  filter(ct < 40) %>%
  filter(obs_time == t) %>% 
  group_by(virus) %>% 
  filter(i %in% sample(unique(i), 100)) %>% 
  ggplot() + 
  geom_jitter(aes(x=confirm_delay,y=ct,col=virus),height=0,width=0.25) +
  geom_smooth(aes(x=confirm_delay,y=ct,col=virus)) +
  scale_y_continuous(trans="reverse")
cts_indiv_sub %>% 
  filter(ct < 40) %>%
  filter(obs_time == t) %>%
  group_by(obs_time,virus) %>% 
  summarize(mean_ct=mean(ct),N=n()) %>% 
  filter(N >= 20) %>%
  ggplot() + geom_line(aes(x=obs_time,y=mean_ct,col=virus))



plot_indiv_simulated_ct_curve(cts_indiv_sub,1:35,100)

N <- 1000000
tmp <- simulate_individual_level_data_symptomatic(N,virus1_inc, times, chain, parTab, max_vl=16,confirm_delays=extraDistr::rdgamma(N,5, 0.8))
