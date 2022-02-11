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

set.seed(123)

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
pars <- c(
  dt=0.25,
  S_ini=1e7,
  
  k_E = 2, ## Number of exposed compartment stages
  k_I = 2, ## Number of infected compartment stages
  sigma = 1/3, ## Latent period
  gamma = 1/7, ## Infectious period
  strain12_immunity = 0.75, ## Cross immunity strain 1 confers against strain 2
  strain21_immunity = 1, ## Cross immunity strain 2 confers against strain 1
  
  beta1=1.5/7, ## Transmission rate in first period
  beta2=1.5/7, ## Transmission rate in second period
  beta3=1.5/7, ## Transmission rate in third period
  beta4=1.5/7, ## Transmission rate in fourth period
  
  tdur1=1000,
  tdur2=1000,
  tdur3=1000,
  tdur4=1000,
  beta_transition_dur=7.25,
  
  strain1_trans=1,
  strain2_trans=4/1.5,
  strain12_trans=0.5,
  strain21_trans=0.5,
  
  wane_rate=1e-05,
  
  importtime1=0,
  importtime2=180,
  seed_size1=100,
  seed_size2=100,
  seed_dur1=7,
  seed_dur2=7
  ) ## Duration of importations

## R0

times <- seq(0, 365*1.5,by=1) ## Run model for 1.5 years

## Change this to change the "sample at day 270" in figure 1 and beyond
samp_time <- 235

seir_dynamics <- run_2strain_seir_simulation(pars, states,times,model_ver="odin",
                                             seir_filename="~/Documents/GitHub/variant_viral_loads/code/odin_files/odinseirstrains.R",
                                             n_repeats=1000)
seir_dynamics$p_inc + theme_bw() + geom_vline(xintercept=samp_time)
virus1_inc <- seir_dynamics$virus1_inc/pars["S_ini"]
virus2_inc <- seir_dynamics$virus2_inc/pars["S_ini"]


virus1_inc_run <- seir_dynamics$virus1_inc_repeats %>% filter(run==1) %>% pull(inc)
virus1_inc_run <- virus1_inc_run/pars["S_ini"]
virus2_inc_run <- seir_dynamics$virus2_inc_repeats %>% filter(run==1) %>% pull(inc)
virus2_inc_run <- virus2_inc_run/pars["S_ini"]

virus_inc <- virus1_inc + virus2_inc

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

p_ct_model <- plot_simulated_ct_curve(vl_pars, ages, 100)
p_ct_model_2 <- plot_simulated_ct_curve_2variants(vl_pars,vl_pars_both,ages,20)

p_ct_compare1 <- p_sim_ct_compare_naive(vl_pars,vl_pars_both,virus1_inc_run,virus2_inc_run, ages,samp_time=samp_time,N=100)
p_ct_compare1
## Using virosolver package, get predicted Ct distribution on each day of the simulation
ct_dist_1 <- calculate_ct_distribution(vl_pars, ages, virus1_inc,times[times >= pars["importtime1"]+7]) %>% mutate(virus="Original variant")
ct_dist_2 <- calculate_ct_distribution(vl_pars, ages, virus2_inc,times[times >= pars["importtime2"]+7]) %>% mutate(virus="New variant, same kinetics")
ct_dist_2_alt <- calculate_ct_distribution(vl_pars_both, ages, virus2_inc,times[times >= pars["importtime2"]+7]) %>% mutate(virus="New variant, different kinetics")
ct_dist_overall <- calculate_ct_distribution(vl_pars, ages, virus_inc,times[times >= pars["importtime1"]+7]) %>% mutate(virus="Overall")
ct_combined_summaries <- bind_rows(ct_dist_1,ct_dist_2,ct_dist_overall,ct_dist_2_alt)
ct_combined_summaries <- ct_combined_summaries %>% left_join(grs) %>% ungroup()

p1 <- plot_medians_and_skew(ct_combined_summaries)


dat_inc_for_plot <- seir_dynamics$seir_res %>% filter(compartment == "obs_inc") %>% 
  mutate(strain=ifelse(strain %in% c("New variant","Reinfection (original->new)"),"New variant","Original variant")) %>%
  group_by(strain, time) %>% summarize(y=mean(y)) %>% select(strain, time, y)
dat_inc_for_plot <- dat_inc_for_plot %>% bind_rows(dat_inc_for_plot %>% group_by(time) %>% summarize(y = sum(y)) %>% mutate(strain="Overall"))
dat_inc_for_plot$strain <- factor(dat_inc_for_plot$strain, levels=c("Overall","Original variant","New variant"))
p_inc <- ggplot(dat_inc_for_plot) + 
  geom_line(aes(x=time,y=y/pars["S_ini"],col=strain))+
  geom_vline(xintercept=c(180),linetype="dashed",col="#D55E00") +
  geom_vline(xintercept=c(0),linetype="dashed",col="#0072B2") +
  variant_color_scale +
  xlab("Time") +
  scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
  scale_y_continuous(breaks=seq(0,0.008,by=0.002))+
  #scale_y_continuous(expand=c(0,0)) +
  ylab("Per capita incidence") +
  theme_overall + 
  theme(legend.position=c(0.8,0.8)) +
  theme_nice_axes+ theme_no_x_axis + labs(tag="A")+ variant_color_scale_fig1 + 
  scale_x_continuous(limits=c(0,425))


p_LHS <- p_inc /
  (p1[[1]] + 
     scale_x_continuous(limits=c(0,425),breaks=seq(0,425,by=50))+ 
     scale_y_continuous(breaks=seq(26,37,by=1)) +
     coord_cartesian(ylim=c(36,26)) +
     labs(tag="C") + 
     geom_vline(xintercept=samp_time,linetype="dotted",col="grey40",size=0.75) + 
     scale_y_continuous(trans="reverse"))
p_RHS <- (p_ct_model_2 + labs(tag="B"))/
  (p_ct_compare1 + labs(tag="D"))

fig1 <- (p_LHS | p_RHS) + plot_layout(widths=c(1.5,1))
ggsave(fig1,filename = "figures/fig1.pdf",height=5,width=8)
ggsave(fig1,filename = "figures/fig1.png",height=5,width=8,dpi=300,units="in")

###############################################
## FIGURE 2
###############################################
gr_tests <- c(0.03,-0.03)
p_aligned <- plot_growth_rate_lineups(ct_combined_summaries %>% mutate(gr=gr_rollmean7) %>% filter(inc > 1e-6, gr > -0.05, gr < 0.05))
p_aligned_median <- p_aligned[[1]] + geom_vline(xintercept=gr_tests,linetype="dotted",col="grey40",size=0.75)
p_ct_samp_gr1 <- p_sim_ct_compare_growth(vl_pars,vl_pars_both,virus1_inc_run,virus2_inc_run, 
                                         ages,ct_combined_summaries,gr_tests[1],N=100,dotsize=1)
p_ct_samp_gr2 <- p_sim_ct_compare_growth(vl_pars,vl_pars_both,virus1_inc_run,virus2_inc_run, 
                                         ages,ct_combined_summaries,gr_tests[2],N=100,dotsize=1)

fig2 <- (p_aligned_median+ theme(legend.position=c(0.85,0.85)) + labs(tag="A")) / 
           (p_ct_samp_gr1+ labs(tag="B"))  / 
           (p_ct_samp_gr2+ labs(tag="C")) 
ggsave(fig2,filename = "figures/fig2.pdf",height=7,width=5)
ggsave(fig2,filename = "figures/fig2.png",height=7,width=5,dpi=300,units="in")


###############################################
## POWER CALCULATION FOR RANDOM CROSS-SECTIONS
###############################################
samp_time <- 235
samp_sizes <- c(25,50,100,250,500)
N_trials <- 100

power_pop_all <- p_sim_ct_compare_power(vl_pars,vl_pars_both,seir_dynamics$virus1_inc_repeats,seir_dynamics$virus2_inc_repeats, 
                                        ages,samp_time=samp_time,trials=N_trials,samp_sizes=samp_sizes)

## Same again but aligned by growth rate
power_pop_gr_up <- p_sim_ct_compare_power(vl_pars,vl_pars_both,seir_dynamics$virus1_inc_repeats,seir_dynamics$virus2_inc_repeats, 
                                          ages,samp_time=samp_time,trials=N_trials,samp_sizes=samp_sizes,
                                   align_gr=TRUE,grs=grs,gr_test=0.03)
power_pop_gr_down <- p_sim_ct_compare_power(vl_pars,vl_pars_both,seir_dynamics$virus1_inc_repeats,seir_dynamics$virus2_inc_repeats, 
                                            ages,samp_time=samp_time,trials=N_trials,samp_sizes=samp_sizes,
                                        align_gr=TRUE,grs=grs,gr_test=-0.03)
ggsave(filename="figures/figS1.png",power_pop_all[[3]],height=8,width=6,units="in",dpi=300)
ggsave(filename="figures/figS2.png",power_pop_gr_up[[3]],height=8,width=6,units="in",dpi=300)
ggsave(filename="figures/figS3.png",power_pop_gr_down[[3]],height=8,width=6,units="in",dpi=300)


###############################################
## Symptomatic reporting population dataset
###############################################
## Using virosolver package, get predicted Ct distribution on each day of the simulation
ct_dist_1_symptomatic <- calculate_ct_distribution(vl_pars, ages, virus1_inc,times[times >= pars["importtime1"]+7],symptom_surveillance = TRUE) %>% mutate(virus="Original variant")
ct_dist_2_symptomatic <- calculate_ct_distribution(vl_pars, ages, virus2_inc,times[times >= pars["importtime2"]+7],symptom_surveillance = TRUE) %>% mutate(virus="New variant, same kinetics")
ct_dist_2_alt_symptomatic <- calculate_ct_distribution(vl_pars_both, ages, virus2_inc,times[times >= pars["importtime2"]+7],symptom_surveillance = TRUE) %>% mutate(virus="New variant, different kinetics")
ct_dist_overall_symptomatic <- calculate_ct_distribution(vl_pars, ages, virus_inc,times[times >= pars["importtime1"]+7],symptom_surveillance = TRUE) %>% mutate(virus="Overall")
ct_combined_summaries_symptomatic <- bind_rows(ct_dist_1_symptomatic,ct_dist_2_symptomatic,ct_dist_overall_symptomatic,ct_dist_2_alt_symptomatic)
ct_combined_summaries_symptomatic <- ct_combined_summaries_symptomatic %>% left_join(grs) %>% ungroup()

p1_symptomatic <- plot_medians_and_skew(ct_combined_summaries_symptomatic %>% filter(inc > 1e-6))

## Plot Ct values over time since infection
p_sympt_ct_kinetics_pop <- plot_simulated_ct_curve_2variants_symptomatic(vl_pars,vl_pars_both, N=1000,xmax=15)
p_sympt_ct_kinetics_pop <- p_sympt_ct_kinetics_pop + labs(tag="A")

p_symptom_compare <- p1_symptomatic[[1]] + 
  coord_cartesian(ylim=c(28.5,21.8)) +
  scale_y_continuous(breaks=seq(22,28,by=1)) +
  geom_vline(xintercept=samp_time,linetype="dotted",size=1,col="grey40")+
  labs(tag="B")

## Plot Ct values over entire epidemic
p_onset_dist <- plot_time_since_infection(samp_time,vl_pars,vl_pars_both,virus1_inc,virus2_inc)  + labs(tag="C")

## Single trial of drawing Ct values and comparing them
p_ct_compare_symp <- p_sim_ct_compare_naive(vl_pars,vl_pars_both,virus1_inc,virus2_inc, ages,samp_time=samp_time,N=100,symptom_surveillance = TRUE) + labs(tag="D")

## Create figure 3
fig3 <- (p_sympt_ct_kinetics_pop | p_symptom_compare )/(p_onset_dist | p_ct_compare_symp)

ggsave(filename="figures/fig3.png",fig3,height=6,width=8,units="in",dpi=300)
ggsave(filename="figures/fig3.pdf",fig3,height=6,width=8)

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

ggsave(filename="figures/figS4.png",power_pop_sympt_lm[[6]],height=7,width=5.5,units="in",dpi=300)
ggsave(filename="figures/figS5.png",power_pop_sympt_lm[[5]],height=7,width=5.5,units="in",dpi=300)
ggsave(filename="figures/figS4.pdf",power_pop_sympt_lm[[6]],height=7,width=5.5)
ggsave(filename="figures/figS5.pdf",power_pop_sympt_lm[[5]],height=7,width=5.5)



