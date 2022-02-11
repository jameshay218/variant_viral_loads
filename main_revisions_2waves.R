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
library(virosolver)
library(ggpubr)
library(lazymcmc)
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

model_pars <- read.csv("pars/partab_seir_model.csv")
vl_pars <- model_pars$values
names(vl_pars) <- model_pars$names

## Change assumed viral kinetics and sampling parameters here
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

###############################################
## 1) SEIR SIMULATION
###############################################
## SEIR parameters
N <- 1e7
pars <- c(
  dt=0.25,
  S_ini=N,
  
  k_E = 2, ## Number of exposed compartment stages
  k_I = 2, ## Number of infected compartment stages
  sigma = 1/3, ## Latent period
  gamma = 1/7, ## Infectious period
  strain12_immunity = 0.75, ## Cross immunity strain 1 confers against strain 2
  strain21_immunity = 1, ## Cross immunity strain 2 confers against strain 1
  
  beta1=1.5/7, ## Transmission rate in first period
  beta2=0.8/7, ## Transmission rate in second period
  beta3=0.8/7, ## Transmission rate in third period
  beta4=0.4/7, ## Transmission rate in third period
  
  
  tdur1=120,
  tdur2=80,
  tdur3=90,
  tdur4=1000,
  beta_transition_dur=7.25,
  
  strain1_trans=1,
  strain2_trans=2,
  strain12_trans=0.5,
  strain21_trans=0.5,
  
  wane_rate=1e-05,
  
  importtime1=0,
  importtime2=180,
  seed_size1=N*1e-6,
  seed_size2=N*1e-6,
  seed_dur1=7,
  seed_dur2=7
) ## Duration of importations


times <- seq(0, 365*1.5,by=1) ## Run model for 1.5 years

seir_dynamics <- run_2strain_seir_simulation(pars, states,times,model_ver="odin",
            seir_filename="~/Documents/GitHub/variant_viral_loads/code/odin_files/odinseirstrains.R",
            n_repeats=1000)

seir_dynamics$p_inc + theme_bw() + geom_vline(xintercept=235)

virus1_inc <- seir_dynamics$virus1_inc/pars["S_ini"]
virus2_inc <- seir_dynamics$virus2_inc/pars["S_ini"]
virus_inc <- virus1_inc + virus2_inc

virus1_inc_run <- seir_dynamics$virus1_inc_repeats %>% filter(run==1) %>% pull(inc)
virus1_inc_run <- virus1_inc_run/pars["S_ini"]
virus2_inc_run <- seir_dynamics$virus2_inc_repeats %>% filter(run==1) %>% pull(inc)
virus2_inc_run <- virus2_inc_run/pars["S_ini"]

## Get all runs
virus1_inc_repeats <- seir_dynamics$virus1_inc_repeats %>% mutate(virus="Original variant") %>% rename(t=time)
virus2_inc_repeats <- seir_dynamics$virus2_inc_repeats%>% mutate(virus="New variant, same kinetics") %>% rename(t=time)
virus2b_inc_repeats <- seir_dynamics$virus2_inc_repeats%>% mutate(virus="New variant, different kinetics") %>% rename(t=time)
virus_inc_repeats <- bind_rows(virus1_inc_repeats,virus2_inc_repeats,virus2b_inc_repeats)


## Pulling incidence curves and getting daily growth rates using repeat data
virus1_inc_repeats_determ <- seir_dynamics$virus1_inc_repeats_determ %>% mutate(virus="Original variant") %>% rename(t=time)
virus2_inc_repeats_determ <- seir_dynamics$virus2_inc_repeats_determ%>% mutate(virus="New variant, same kinetics") %>% rename(t=time)
virus2b_inc_repeats_determ <- seir_dynamics$virus2_inc_repeats_determ%>% mutate(virus="New variant, different kinetics") %>% rename(t=time)
virus_inc_repeats_determ <- bind_rows(virus1_inc_repeats_determ,virus2_inc_repeats_determ,virus2b_inc_repeats_determ)

gr_repeats <- virus_inc_repeats_determ %>% group_by(run,virus) %>% 
  mutate(gr = log(lead(inc,1)/inc)) %>%
  mutate(gr=ifelse(is.infinite(gr) | is.nan(gr),0,gr)) %>%
  mutate(gr_rollmean7=zoo::rollmean(gr,k=7,fill=NA))

ggplot(gr_repeats) + geom_line(aes(x=t,y=gr_rollmean7,group=interaction(run,virus),col=virus),size=0.1) + coord_cartesian(ylim=c(-0.25,0.25))

gr1 <- tibble(t=times,gr=c(0,log(virus1_inc[2:length(virus1_inc)]/virus1_inc[1:(length(virus1_inc)-1)])),inc=virus1_inc,virus="Original variant")
gr2 <- tibble(t=times,gr=c(0,log(virus2_inc[2:length(virus2_inc)]/virus2_inc[1:(length(virus2_inc)-1)])),inc=virus2_inc,virus="New variant, same kinetics")
gr2_alt <- tibble(t=times,gr=c(0,log(virus2_inc[2:length(virus2_inc)]/virus2_inc[1:(length(virus2_inc)-1)])),inc=virus2_inc,virus="New variant, different kinetics")
gr_overall <- tibble(t=times,gr=c(0,log(virus_inc[2:length(virus_inc)]/virus_inc[1:(length(virus_inc)-1)])),inc=virus_inc,virus="Overall")
grs <- bind_rows(gr1, gr2,gr2_alt,gr_overall)

grs <- grs %>% dplyr::arrange(virus,t) %>% group_by(virus) %>% 
  mutate(gr=ifelse(is.infinite(gr),0,gr)) %>%
  mutate(gr_rollmean7=zoo::rollmean(gr,k=7,fill=NA))

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
## Change this to change the "sample at day 270" in figure 1 and beyond
samp_time <- 235

p_ct_model <- plot_simulated_ct_curve(vl_pars, ages, 100)
p_ct_model_2 <- plot_simulated_ct_curve_2variants(vl_pars,vl_pars_both,ages,20)

p_ct_compare1 <- p_sim_ct_compare_naive(vl_pars,vl_pars_both,virus1_inc,virus2_inc, ages,samp_time=samp_time,N=100)

## Using virosolver package, get predicted Ct distribution on each day of the simulation
ct_dist_1 <- calculate_ct_distribution(vl_pars, ages, virus1_inc,times[times >= pars["importtime1"]+7]) %>% mutate(virus="Original variant")
ct_dist_2 <- calculate_ct_distribution(vl_pars, ages, virus2_inc,times[times >= pars["importtime2"]+7]) %>% mutate(virus="New variant, same kinetics")
ct_dist_2_alt <- calculate_ct_distribution(vl_pars_both, ages, virus2_inc,times[times >= pars["importtime2"]+7]) %>% mutate(virus="New variant, different kinetics")
ct_dist_overall <- calculate_ct_distribution(vl_pars, ages, virus_inc,times[times >= pars["importtime1"]+7]) %>% mutate(virus="Overall")
ct_combined_summaries <- bind_rows(ct_dist_1,ct_dist_2,ct_dist_overall,ct_dist_2_alt)
ct_combined_summaries <- ct_combined_summaries %>% left_join(grs) %>% ungroup()

p1 <- plot_medians_and_skew(ct_combined_summaries)

beta_dat <- seir_dynamics$seir_res %>% filter(compartment == "beta")
use_grs <- grs %>% filter(inc > 1e-5) %>%
  filter(virus %in% c("New variant, same kinetics","Original variant","Overall")) %>%
  mutate(virus = ifelse(virus == "New variant, same kinetics", "New variant",virus))
use_grs$virus<- factor(use_grs$virus, levels=c("Overall","New variant","Original variant"))
p_gr_new <- ggplot(use_grs) +
  geom_line(aes(x=t,y=gr_rollmean7,col=virus)) + 
  coord_cartesian(ylim=c(-0.07,0.07)) +
  variant_color_scale_stochastic_alt + geom_vline(xintercept=samp_time,linetype="dotted",col="grey40",size=0.75)+
  #geom_vline(xintercept=c(75,250),linetype="dotted",col="red",size=0.75) +
  #geom_vline(xintercept=c(150,350),linetype="dotted",col="darkgreen",size=0.75) +
  theme_overall + theme_nice_axes +
  ylab("7-day rolling average growth rate") +
  theme(legend.position=c(0.9,0.8)) +
  scale_x_continuous(limits=c(0,425)) + 
  labs(tag="C")
p_gr_new

p_inc_new <- seir_dynamics$p_inc + 
  geom_line(data=beta_dat,aes(x=time,y=y/400),col="black",size=0.75,alpha=0.5) +
  labs(tag="A") + 
  scale_x_continuous(limits=c(0,425)) + 
  scale_y_continuous(expand=c(1e-4,1e-4),name="Incidence",
                     breaks=seq(0,0.003,by=0.0005),
                     sec.axis=sec_axis(trans=~.*400,name="Daily contact rate (grey)")) +
  variant_color_scale_stochastic + variant_linetype_scale_stochastic + 
  theme(legend.position=c(0.9,0.8)) +
  geom_vline(xintercept=samp_time,linetype="dotted",col="grey40",size=0.75)

p_LHS <- (p_inc_new)/
  (p1[[1]]+ scale_x_continuous(limits=c(0,425)) + labs(tag="B") + 
     scale_y_continuous(trans="reverse") + theme_overall + theme_nice_axes + theme_no_x_axis +
     theme(legend.position=c(0.9,0.8)) + geom_vline(xintercept=samp_time,linetype="dotted",col="grey40",size=0.75))/
  p_gr_new

fig1 <- p_LHS
ggsave(fig1,filename = "figures_revisions/new_figS1.pdf",height=8,width=6)
ggsave(fig1,filename = "figures_revisions/new_figS1.png",height=8,width=6,dpi=300,units="in")

###############################################
## POWER CALCULATION FOR RANDOM CROSS-SECTIONS
###############################################
samp_time <- 235
samp_sizes <- c(25,50,100,250,500)
N_trials <- 100

power_pop_all <- p_sim_ct_compare_power(vl_pars,vl_pars_both,seir_dynamics$virus1_inc_repeats,
                                        seir_dynamics$virus2_inc_repeats, ages,samp_time=samp_time,
                                        trials=N_trials,samp_sizes=samp_sizes)

## Same again but aligned by growth rate
## Because of the non-smooth growth rate transitions, we are going to home in on periods of
## stable growth
power_pop_gr_up <- p_sim_ct_compare_power(vl_pars,vl_pars_both,seir_dynamics$virus1_inc_repeats,
                                          seir_dynamics$virus2_inc_repeats, ages,samp_time=samp_time,trials=N_trials,samp_sizes=samp_sizes,
                                          align_gr=TRUE,grs=grs,gr_test=0.06)

grs_use_down <- grs %>% filter((virus == "Original variant" & t > 150 & t < 200) | 
                                 (virus %in% c("New variant, different kinetics","New variant, same kinetics") & t > 325 & t < 400))
power_pop_gr_down <- p_sim_ct_compare_power(vl_pars,vl_pars_both,seir_dynamics$virus1_inc_repeats,
                                            seir_dynamics$virus2_inc_repeats, ages,samp_time=samp_time,trials=N_trials,samp_sizes=samp_sizes,
                                            align_gr=TRUE,grs=grs_use_down,gr_test=-0.025)

ggsave(filename="figures_revisions/new_figS2.png",power_pop_all[[3]],height=8,width=6,units="in",dpi=300)
ggsave(filename="figures_revisions/new_figS3.png",power_pop_gr_up[[3]],height=8,width=6,units="in",dpi=300)
ggsave(filename="figures_revisions/new_figS4.png",power_pop_gr_down[[3]],height=8,width=6,units="in",dpi=300)

ggsave(filename="figures_revisions/new_figS2.pdf",power_pop_all[[3]],height=8,width=6,units="in")
ggsave(filename="figures_revisions/new_figS3.pdf",power_pop_gr_up[[3]],height=8,width=6,units="in")
ggsave(filename="figures_revisions/new_figS4.pdf",power_pop_gr_down[[3]],height=8,width=6,units="in")


###############################################
## Symptomatic reporting population dataset
###############################################
## Using virosolver package, get predicted Ct distribution on each day of the simulation
ct_dist_1_symptomatic <- calculate_ct_distribution(vl_pars, ages, virus1_inc,
                                                   times[times >= pars["importtime1"]+35],
                                                   symptom_surveillance = TRUE) %>% mutate(virus="Original variant")

ct_dist_2_symptomatic <- calculate_ct_distribution(vl_pars, ages, virus2_inc,times[times >= pars["importtime2"]+35],symptom_surveillance = TRUE) %>% mutate(virus="New variant, same kinetics")
ct_dist_2_alt_symptomatic <- calculate_ct_distribution(vl_pars_both, ages, virus2_inc,times[times >= pars["importtime2"]+35],symptom_surveillance = TRUE) %>% mutate(virus="New variant, different kinetics")
ct_dist_overall_symptomatic <- calculate_ct_distribution(vl_pars, ages, virus_inc,times[times >= pars["importtime1"]+35],symptom_surveillance = TRUE) %>% mutate(virus="Overall")
ct_combined_summaries_symptomatic <- bind_rows(ct_dist_1_symptomatic,ct_dist_2_symptomatic,ct_dist_overall_symptomatic,ct_dist_2_alt_symptomatic)
ct_combined_summaries_symptomatic <- ct_combined_summaries_symptomatic %>% left_join(grs) %>% ungroup()

p1_symptomatic <- plot_medians_and_skew(ct_combined_summaries_symptomatic)

p_symptom_compare <- p1_symptomatic[[1]] + 
  coord_cartesian(ylim=c(30.5,20.8)) +
  scale_y_continuous(breaks=seq(21,30,by=1)) +
  scale_x_continuous(limits=c(0,300)) + labs(tag="B") + 
  theme(legend.position=c(0.18,0.78)) + geom_vline(xintercept=samp_time,linetype="dotted",col="grey40",size=0.75)
  geom_vline(xintercept=samp_time,linetype="dotted",size=1,col="grey40")

## Create figure 3
figS5 <- p_inc_new/p_symptom_compare

ggsave(filename="figures_revisions/new_figS5.png",figS5,height=5,width=6,units="in",dpi=300)
ggsave(filename="figures_revisions/new_figS5.pdf",figS5,height=5,width=6)

## Similar simulation, but using a re-weighted incidence curve around the sample time to ensure that we have plenty of Ct values
## to resample from. Using the above version of the simulation generates too few Ct values to be reliable.
ct_dist_symptomatic_window <- simulate_popn_cts_symptomatic_repeats(virus1_inc_repeats, 
                                                                    virus2_inc_repeats,
                                                                    samp_time, 30, 10,
                                                            vl_pars, vl_pars_both,
                                                            pars["S_ini"], 
                                                            confirm_delay_par1_v1=5, confirm_delay_par2_v1 = 1,
                                                            confirm_delay_par1_v2=5, confirm_delay_par2_v2 = 1)

ct_dist_symptomatic_window_dat <- ct_dist_symptomatic_window$ct_dat_sympt
ct_dist_symptomatic_window_dat <- ct_dist_symptomatic_window_dat %>% mutate(days_since_onset = sampled_time - onset_time)
ct_dist_symptomatic_window_dat <- ct_dist_symptomatic_window_dat%>% dplyr::select(-ct) %>% rename(ct=ct_obs)

## Power calculation comparing the use of a linear regression model and Wilcoxon rank sum test
power_pop_sympt_lm <- p_sim_ct_compare_power_symp_regression(ct_dist_symptomatic_window_dat %>% filter(ct < 40), 
                            samp_time=samp_time,trials=N_trials,samp_size=samp_sizes,
                            true_peak_diff=vl_pars_both["viral_peak"]-vl_pars["viral_peak"],samp_window=7)

ggsave(filename="figures_revisions/new_figS6.png",power_pop_sympt_lm[[6]] + theme(legend.position=c(0.8,0.4)),height=7,width=5.5,units="in",dpi=300)
ggsave(filename="figures_revisions/new_figS7.png",power_pop_sympt_lm[[5]],height=7,width=5.5,units="in",dpi=300)
ggsave(filename="figures_revisions/new_figS6.pdf",power_pop_sympt_lm[[6]]+ theme(legend.position=c(0.8,0.4)),height=7,width=5.5)
ggsave(filename="figures_revisions/new_figS7.pdf",power_pop_sympt_lm[[5]],height=7,width=5.5)



