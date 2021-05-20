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
#devtools::load_all("~/Documents/GitHub/virosolver")
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
p_ct_model <- plot_simulated_ct_curve(vl_pars, ages, 100)
p_ct_model_2 <- plot_simulated_ct_curve_2variants(vl_pars,vl_pars_both,ages,20)

p_ct_compare1 <- p_sim_ct_compare_naive(vl_pars,vl_pars_both,virus1_inc,virus2_inc, ages,samp_time=270,N=100)

## Using virosolver package, get predicted Ct distribution on each day of the simulation
ct_dist_1 <- calculate_ct_distribution(vl_pars, ages, virus1_inc,times[times >= pars["importtime1"]+25]) %>% mutate(virus="Original variant")
ct_dist_2 <- calculate_ct_distribution(vl_pars, ages, virus2_inc,times[times >= pars["importtime2"]+25]) %>% mutate(virus="New variant, same kinetics")
ct_dist_2_alt <- calculate_ct_distribution(vl_pars_both, ages, virus2_inc,times[times >= pars["importtime2"]+25]) %>% mutate(virus="New variant, different kinetics")
ct_dist_overall <- calculate_ct_distribution(vl_pars, ages, virus_inc,times[times >= pars["importtime1"]+25]) %>% mutate(virus="Overall")
ct_combined_summaries <- bind_rows(ct_dist_1,ct_dist_2,ct_dist_overall,ct_dist_2_alt)
ct_combined_summaries <- ct_combined_summaries %>% left_join(grs) %>% ungroup()

p1 <- plot_medians_and_skew(ct_combined_summaries)

p_LHS <- (seir_dynamics$p_inc + labs(tag="A"))/
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
           (p_ct_samp_gr1+ labs(tag="A"))  / 
           (p_ct_samp_gr2+ labs(tag="A")) 
ggsave(fig2,filename = "figures/fig2.pdf",height=7,width=5)
ggsave(fig2,filename = "figures/fig2.png",height=7,width=5,dpi=300,units="in")

