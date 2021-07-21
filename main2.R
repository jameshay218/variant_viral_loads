###############################################
## CT VALUES UNDER COMPETING STRAIN DYNAMICS
## James Hay
## jhay@hsph.harvard.edu
## May 19, 2021
###############################################

###############################################
## HEADERS
###############################################
library(tidyverse)
library(patchwork)
library(extraDistr)
#library(virosolver)
library(ggpubr)
library(lazymcmc)
library(rethinking)
devtools::load_all("~/Documents/GitHub/virosolver")
## Where to perform the simulations
HOME_WD <- "~"
HOME_WD <- "~/Documents/GitHub/variant_viral_loads/"
setwd(HOME_WD)

## Load functions for line list simulation
source("code/linelist_sim_funcs.R")
source("code/plotting.R")
source("code/seir_funcs.R")
source("code/analysis_funcs.R")
source("code/invasion_rates_KISSLER2020.R")
source("code/simulate_symptomatic_population.R")
source('~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/functions/simulation_functions.R')
source("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/functions/model_funcs_multivariate_hinge.R")


model_pars <- read.csv("pars/partab_seir_model.csv")
vl_pars <- model_pars$values
names(vl_pars) <- model_pars$names

vl_pars["incu_par1"] <- 1.621
vl_pars["incu_par2"] <- 0.418

vl_pars["sampling_par1"] <- 5
vl_pars["sampling_par2"] <- 1

vl_pars["max_incu_period"] <- 25
vl_pars["max_sampling_delay"] <- 25

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

## MCMC chains for individual-level kinetics
parTab <- read.csv("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/pars/partab_multivariate_hinge.csv",stringsAsFactors=FALSE)
chains <- load_mcmc_chains(paste0("~/Documents/GitHub/covid19-group-tests/code/viral_kinetics/chains/chains_swab"), parTab, 
                           FALSE, 100, burnin=1000000,multi=TRUE)
chain <- as.data.frame(chains[[2]])
chain$sampno <- 1:nrow(chain)

chain_vl2 <- chain
chain_vl2$viral_peak_mean <- chain_vl2$viral_peak_mean + 1
chain_vl2$wane_mean <- chain_vl2$wane_mean + 5


###############################################
## 1) SEIR SIMULATION
###############################################
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

gr1 <- tibble(t=times,gr=c(0,log(virus1_inc[2:length(virus1_inc)]/virus1_inc[1:(length(virus1_inc)-1)])),inc=virus1_inc,virus="Original variant")
gr2 <- tibble(t=times,gr=c(0,log(virus2_inc[2:length(virus2_inc)]/virus2_inc[1:(length(virus2_inc)-1)])),inc=virus_inc,virus="New variant, same kinetics")
gr2_alt <- tibble(t=times,gr=c(0,log(virus2_inc[2:length(virus2_inc)]/virus2_inc[1:(length(virus2_inc)-1)])),inc=virus_inc,virus="New variant, different kinetics")
gr_overall <- tibble(t=times,gr=c(0,log(virus_inc[2:length(virus_inc)]/virus_inc[1:(length(virus_inc)-1)])),inc=virus_inc,virus="Overall")
grs <- bind_rows(gr1, gr2,gr2_alt,gr_overall)

## Assume individuals can stay detectable for up to 35 days
lastday <- 35
ages <- 1:lastday

## Calculate viral load model and proportion detectable over time based on Hay & Kennedy-Shaffer et al. 2021
viral_loads <- viral_load_func(vl_pars, ages) ## Get modal Ct
detectable_props <- prop_detectable(ages,vl_pars,viral_loads) ## Get proportion detectable over days since infection

## Get infection age distributions over time WRT virus 1, virus 2 and overall
age_dist_v1 <- calculate_infection_age_distribution(seir_dynamics$virus1_inc,detectable_props,ages)
age_dist_v2 <- calculate_infection_age_distribution(seir_dynamics$virus2_inc,detectable_props,ages)
age_dist_comb <- calculate_infection_age_distribution(seir_dynamics$virus1_inc+seir_dynamics$virus2_inc,detectable_props,ages)

###############################################
## FIGURE 1
###############################################
samp_time <- 270

p_ct_model <- plot_simulated_ct_curve(vl_pars, ages, 100)
p_ct_model_2 <- plot_simulated_ct_curve_2variants(vl_pars,vl_pars_both,ages,20)

p_ct_compare1 <- p_sim_ct_compare_naive(vl_pars,vl_pars_both,virus1_inc,virus2_inc, ages,samp_time=samp_time,N=100)

## Using virosolver package, get predicted Ct distribution on each day of the simulation
ct_dist_1 <- calculate_ct_distribution(vl_pars, ages, virus1_inc,times[times >= pars["importtime1"]+35]) %>% mutate(virus="Original variant")
ct_dist_2 <- calculate_ct_distribution(vl_pars, ages, virus2_inc,times[times >= pars["importtime2"]+35]) %>% mutate(virus="New variant, same kinetics")
ct_dist_2_alt <- calculate_ct_distribution(vl_pars_both, ages, virus2_inc,times[times >= pars["importtime2"]+35]) %>% mutate(virus="New variant, different kinetics")
ct_dist_overall <- calculate_ct_distribution(vl_pars, ages, virus_inc,times[times >= pars["importtime1"]+35]) %>% mutate(virus="Overall")
ct_combined_summaries <- bind_rows(ct_dist_1,ct_dist_2,ct_dist_overall,ct_dist_2_alt)
ct_combined_summaries <- ct_combined_summaries %>% left_join(grs) %>% ungroup()

p1 <- plot_medians_and_skew(ct_combined_summaries)

p_LHS <- (seir_dynamics$p_inc + labs(tag="A")+ variant_color_scale_fig1)/
  (p1[[1]] + labs(tag="C") + geom_vline(xintercept=samp_time,linetype="dotted",col="grey40",size=0.75) + 
     scale_y_continuous(trans="reverse"))
p_RHS <- (p_ct_model_2 + labs(tag="B"))/
  (p_ct_compare1 + labs(tag="D"))

fig1 <- (p_LHS | p_RHS) + plot_layout(widths=c(1.5,1))
ggsave(fig1,filename = "figures/fig1.pdf",height=5,width=8)
ggsave(fig1,filename = "figures/fig1.png",height=5,width=8,dpi=300,units="in")

###############################################
## FIGURE 2
###############################################
gr_tests <- c(0.03,-0.02)
p_aligned <- plot_growth_rate_lineups(ct_combined_summaries)
p_aligned_median <- p_aligned[[1]] + geom_vline(xintercept=gr_tests,linetype="dotted",col="grey40",size=0.75)
set.seed(2)
p_ct_samp_gr1 <- p_sim_ct_compare_growth(vl_pars,vl_pars_both,virus1_inc,virus2_inc, ages,combined_summaries,gr_tests[1],N=100,dotsize=1)
set.seed(7)
p_ct_samp_gr2 <- p_sim_ct_compare_growth(vl_pars,vl_pars_both,virus1_inc,virus2_inc, ages,combined_summaries,gr_tests[2],N=100,dotsize=1)

fig2 <- (p_aligned_median + labs(tag="A")) / 
           (p_ct_samp_gr1+ labs(tag="B"))  / 
           (p_ct_samp_gr2+ labs(tag="C")) 
ggsave(fig2,filename = "figures/fig2.pdf",height=7,width=5)
ggsave(fig2,filename = "figures/fig2.png",height=7,width=5,dpi=300,units="in")

###############################################
## FIGURE 1 -- INDIVIDUAL-LEVEL KINETICS
###############################################
samp_size <- 100000
ct_dist_1_indiv <- simulate_individual_level_data(samp_size, virus1_inc, times, chain, parTab, max_vl=16) %>% mutate(virus="Original variant")
ct_dist_2_indiv <- simulate_individual_level_data(samp_size, virus2_inc, times, chain, parTab, max_vl=16) %>% mutate(virus="New variant, same kinetics") %>% mutate(i = i + max(ct_dist_1_indiv$i))
ct_dist_2_indiv_alt <- simulate_individual_level_data(samp_size, virus2_inc, times, chain_vl2, parTab, max_vl=16) %>% mutate(virus="New variant, different kinetics") %>% mutate(i = i + max(ct_dist_2_indiv$i))
cts_indiv_comb <- bind_rows(ct_dist_1_indiv,ct_dist_2_indiv,ct_dist_2_indiv_alt) %>% mutate(virus=factor(virus,levels=variant_levels))
cts_indiv_comb <- cts_indiv_comb %>% filter(ct < 40)
cts_indiv_summaries <- calculate_ct_dist_indiv_summaries(cts_indiv_comb) %>% left_join(grs) %>% ungroup()
p2 <- plot_medians_and_skew(cts_indiv_summaries %>% filter(N > 25))
p_vl_curves_indiv <- plot_indiv_simulated_ct_curve(cts_indiv_comb%>% filter(days_since_infection >= 0,days_since_infection < 35),1:35,100)
p_grs_indiv <- plot_growth_rate_lineups(cts_indiv_summaries%>%filter(N > 50))
fig1_alt <- (p_vl_curves_indiv+labs(tag="A"))/
  (p2[[1]]+labs(tag="B"))/
  (p_grs_indiv[[1]]+labs(tag="C"))
ggsave(fig1_alt,filename = "figures/fig1_alt.pdf",height=7,width=5)
ggsave(fig1_alt,filename = "figures/fig1_alt.png",height=7,width=5,dpi=300,units="in")


###############################################
## POWER CALCULATION FOR RANDOM CROSS-SECTIONS
###############################################
samp_time <- 270
samp_sizes <- c(25,50,100,250,500)
N_trials <- 1000

power_pop_all <- p_sim_ct_compare_power(vl_pars,vl_pars_both,virus1_inc,virus2_inc, ages,samp_time=samp_time,trials=N_trials,samp_sizes=samp_sizes)

## Same again but aligned by growth rate
power_pop_gr_up <- p_sim_ct_compare_power(vl_pars,vl_pars_both,virus1_inc,virus2_inc, ages,samp_time=samp_time,trials=N_trials,samp_sizes=samp_sizes,
                                   align_gr=TRUE,grs=grs,gr_test=0.03)
power_pop_gr_down <- p_sim_ct_compare_power(vl_pars,vl_pars_both,virus1_inc,virus2_inc, ages,samp_time=samp_time,trials=N_trials,samp_sizes=samp_sizes,
                                        align_gr=TRUE,grs=grs,gr_test=-0.02)
ggsave(filename="figures/figS1.png",power_pop_all[[3]],height=8,width=6,units="in",dpi=300)
ggsave(filename="figures/figS2.png",power_pop_gr_up[[3]],height=8,width=6,units="in",dpi=300)
ggsave(filename="figures/figS3.png",power_pop_gr_down[[3]],height=8,width=6,units="in",dpi=300)


peak_ct_vl1 <- 40 - log2(10)*(mean(chain$viral_peak_mean)-2)
peak_ct_vl2 <- 40 - log2(10)*(mean(chain_vl2$viral_peak_mean)-2)
true_diff_indiv_ct <- peak_ct_vl2 -  peak_ct_vl1
## Same again but with individual-level kinetics
power_indiv_all <- p_sim_ct_indiv_compare_power(cts_indiv_comb,270,trials=N_trials,samp_sizes=samp_sizes,alpha=0.05,true_peak_diff = true_diff_indiv_ct)
power_indiv_gr_up <- p_sim_ct_indiv_compare_power(cts_indiv_comb,270,trials=N_trials,samp_sizes=samp_sizes,alpha=0.05,
                                                  align_gr=TRUE,grs=grs,gr_test=0.03,true_peak_diff = true_diff_indiv_ct)
power_indiv_gr_down <- p_sim_ct_indiv_compare_power(cts_indiv_comb,270,trials=N_trials,samp_sizes=samp_sizes,alpha=0.05,
                                                    align_gr=TRUE,grs=grs,gr_test=-0.02,true_peak_diff = true_diff_indiv_ct)

ggsave(filename="figures/figS1_indiv.png",power_indiv_all[[3]],height=7,width=5.5,units="in",dpi=300)
ggsave(filename="figures/figS2_indiv_up.png",power_indiv_gr_up[[3]],height=7,width=5.5,units="in",dpi=300)
ggsave(filename="figures/figS3_indiv_down.png",power_indiv_gr_down[[3]],height=7,width=5.5,units="in",dpi=300)

###############################################
## Symptomatic reporting population dataset
###############################################
## Using virosolver package, get predicted Ct distribution on each day of the simulation
ct_dist_1_symptomatic <- calculate_ct_distribution(vl_pars, ages, virus1_inc,times[times >= pars["importtime1"]+35],symptom_surveillance = TRUE) %>% mutate(virus="Original variant")
ct_dist_2_symptomatic <- calculate_ct_distribution(vl_pars, ages, virus2_inc,times[times >= pars["importtime2"]+35],symptom_surveillance = TRUE) %>% mutate(virus="New variant, same kinetics")
ct_dist_2_alt_symptomatic <- calculate_ct_distribution(vl_pars_both, ages, virus2_inc,times[times >= pars["importtime2"]+35],symptom_surveillance = TRUE) %>% mutate(virus="New variant, different kinetics")
ct_dist_overall_symptomatic <- calculate_ct_distribution(vl_pars, ages, virus_inc,times[times >= pars["importtime1"]+35],symptom_surveillance = TRUE) %>% mutate(virus="Overall")
ct_combined_summaries_symptomatic <- bind_rows(ct_dist_1_symptomatic,ct_dist_2_symptomatic,ct_dist_overall_symptomatic,ct_dist_2_alt_symptomatic)
ct_combined_summaries_symptomatic <- ct_combined_summaries_symptomatic %>% left_join(grs) %>% ungroup()

p1_symptomatic <- plot_medians_and_skew(ct_combined_summaries_symptomatic)

## Plot Ct values over time since infection
p_sympt_ct_kinetics_pop <- plot_simulated_ct_curve_2variants_symptomatic(vl_pars,vl_pars_both, N=1000,xmax=20)
p_sympt_ct_kinetics_pop <- p_sympt_ct_kinetics_pop + labs(tag="A")

p_symptom_compare <- p1_symptomatic[[1]] + 
  coord_cartesian(ylim=c(28.5,21.8)) +
  scale_y_continuous(breaks=seq(22,28,by=1)) +
  geom_vline(xintercept=samp_time,linetype="dotted",size=1,col="grey40")+
  labs(tag="B")

## Plot Ct values over entire epidemic
p_onset_dist <- plot_time_since_infection(270,vl_pars,vl_pars_both,virus1_inc,virus2_inc)  + labs(tag="C")

## Single trial of drawing Ct values and comparing them
p_ct_compare_symp <- p_sim_ct_compare_naive(vl_pars,vl_pars_both,virus1_inc,virus2_inc, ages,samp_time=samp_time,N=100,symptom_surveillance = TRUE) + labs(tag="D")

## Create figure 3
fig3 <- (p_sympt_ct_kinetics_pop | p_symptom_compare )/(p_onset_dist | p_ct_compare_symp)

ggsave(filename="figures/fig3.png",fig3,height=6,width=8,units="in",dpi=300)
ggsave(filename="figures/fig3.pdf",fig3,height=6,width=8)

## Similar simulation, but using a re-weighted incidence curve around the sample time to ensure that we have plenty of Ct values
## to resample from. Using the above version of the simulation generates too few Ct values to be reliable.
ct_dist_symptomatic_window <- simulate_popn_cts_symptomatic(virus1_inc[(samp_time-30):(samp_time+10)]/sum(virus1_inc[(samp_time-30):(samp_time+10)]), 
                                                     virus2_inc[(samp_time-30):(samp_time+10)]/sum(virus2_inc[(samp_time-30):(samp_time+10)]), 
                                                     vl_pars, vl_pars_both,
                                                     1000000, times[(samp_time-30):(samp_time+10)],
                                                     confirm_delay_par1_v1=7, confirm_delay_par2_v1=0.9,
                                                     confirm_delay_par1_v2=5, confirm_delay_par2_v2 = 1)
ct_dist_symptomatic_window_dat <- ct_dist_symptomatic_window$ct_dat_sympt
ct_dist_symptomatic_window_dat <- ct_dist_symptomatic_window_dat %>% mutate(days_since_onset = sampled_time - onset_time)
ct_dist_symptomatic_window_dat <- ct_dist_symptomatic_window_dat%>% dplyr::select(-ct) %>% rename(ct=ct_obs)

## Power calculation comparing the use of a linear regression model and Wilcoxon rank sum test
power_pop_sympt_lm <- p_sim_ct_compare_power_symp_regression(ct_dist_symptomatic_window_dat %>% filter(ct < 40), 
                            samp_time=samp_time,trials=N_trials,samp_size=samp_sizes,
                            true_peak_diff=vl_pars_both["viral_peak"]-vl_pars["viral_peak"],samp_window=7)

ggsave(filename="figures/figS4_diff.png",power_pop_sympt_lm[[6]],height=7,width=5.5,units="in",dpi=300)
ggsave(filename="figures/figS5_diff.png",power_pop_sympt_lm[[5]],height=7,width=5.5,units="in",dpi=300)
ggsave(filename="figures/figS4_diff.pdf",power_pop_sympt_lm[[6]],height=7,width=5.5)
ggsave(filename="figures/figS5_diff.pdf",power_pop_sympt_lm[[5]],height=7,width=5.5)



