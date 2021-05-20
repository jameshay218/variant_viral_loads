###############################################
## CT VALUES UNDER COMPETING STRAIN DYNAMICS
## James Hay
## jhay@hsph.harvard.edu
## May 3, 2021
###############################################
## This script demonstrates how Ct values observed in a population with two, competing circulating strains is expected to vary over time.
## The simulation takes part as:
##    1) Simulate a two-strain SEIR model, where one strain out-competes the other through a combination of immune escape and increased transmissibility
##    2) The distributions of Ct values and infection ages are calculated based on the convolution of the incidence curve and assumed viral kinetics
##    3) Stochastic Ct distributions are simulated based on cross-sectional surveillance, where individual Ct values are simulated
##    4) Stochastic Ct distributions are simulated based on symptom-based surveillance, where individuals are tested if they become symptomatic
##    5) Stochastic Ct distributions are simulated based on a mixture of symptom-based and completely random surveillance
##    6) 2-5 are re-simulated assuming that the new variant has different viral kinetics
##    7) 2-5 are re-simulated assuming that the new variant has only increased transmissibility
##    8) 2-5 are re-simulated assuming that the new variant has only immune escape

###############################################
## HEADERS
###############################################
library(deSolve)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(lazymcmc)
library(extraDistr)
library(virosolver)

## Where to perform the simulations
HOME_WD <- "~"
HOME_WD <- "~/Documents/GitHub/variant_viral_loads/"
setwd(HOME_WD)

## Load functions for line list simulation
source("code/linelist_sim_funcs.R")
source("code/plotting.R")
source("code/invasion_rates_KISSLER2020.R")

set.seed(1)
###############################################

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
## Simulate SEIR model
seir_sim <- as_tibble(as.data.frame(lsoda(states,times,seir_model_2strains,pars)))

## Sense check compartments
seir_sim %>% 
  pivot_longer(-time) %>% 
  filter(!(name %in% c("inc1","inc2"))) %>% 
  ggplot() + geom_line(aes(x=time,y=value,col=name))


## EXTRACT PER CAPITA INCIDENCE FOR EACH VIRUS
virus1_inc <- c(0,diff(seir_sim$inc1))
virus2_inc <- c(0,diff(seir_sim$inc2))

p_inc <- ggplot(seir_sim) + 
  geom_line(aes(x=time,y=c(0,diff(inc1)),col="Original variant"))+
  geom_line(aes(x=time,y=c(0,diff(inc2)),col="New variant")) +
  geom_line(aes(x=time,y=(c(0,diff(inc1)) + c(0,diff(inc2))),col="Overall")) +
  geom_vline(xintercept=c(180),linetype="dashed",col="red") +
  geom_vline(xintercept=c(0),linetype="dashed",col="blue") +
  scale_color_manual(name="Variant",values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  xlab("Time") +
  scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
  #scale_y_continuous(expand=c(0,0)) +
  ylab("Per capita incidence") +
  theme_overall + 
  theme(legend.position=c(0.8,0.8)) +
  theme_nice_axes+ theme_no_x_axis  +
  labs(tag="A")


###############################################
## 2) DETERMINISTIC CT VALUE DISTRIBUTIONS
###############################################
## Viral load model parameters
model_pars <- read.csv("pars/partab_seir_model.csv")
vl_pars <- model_pars$values
names(vl_pars) <- model_pars$names

## Assume individuals can stay detectable for up to 35 days
lastday <- 35
ages <- 1:lastday
viral_loads <- viral_load_func(vl_pars, 1:lastday) ## Get modal Ct
detectable_props <- prop_detectable(1:lastday,vl_pars,viral_loads) ## Get proportion detectable over days since infection

## Find distribution of ages since infection among PCR positive individuals as convolution of incidence curve and proportion detectable
age_res1 <- matrix(nrow=length(virus1_inc),ncol=lastday)
age_res2 <- matrix(nrow=length(virus1_inc),ncol=lastday)
age_res_overall <- matrix(nrow=length(virus1_inc),ncol=lastday)
inc_overall <- virus1_inc + virus2_inc
for (i in 2:length(virus1_inc)) {
  past_inc1 <- virus1_inc[(i-1):(max(i-lastday,1))]
  past_inc2 <- virus2_inc[(i-1):(max(i-lastday,1))]
  past_inc_overall <- inc_overall[(i-1):(max(i-lastday,1))]
  days <- 1:length(past_inc1)
  age_res1[i,days] <- past_inc1*detectable_props[days]
  age_res2[i,days] <- past_inc2*detectable_props[days]
  age_res_overall[i,days] <- past_inc_overall*detectable_props[days]
}

## Using virosolver package, get predicted Ct distribution on each day of the simulation
dat1 <- pred_dist_wrapper(seq(0,40,by=0.01),obs_times = times[times >= pars["importtime1"]+25],ages,vl_pars,virus1_inc)
dat2 <- pred_dist_wrapper(seq(0,40,by=0.01),obs_times = times[times >= pars["importtime2"]+25],ages,vl_pars,virus2_inc)
dat_overall <- pred_dist_wrapper(seq(0,40,by=0.01),obs_times = times[times >= pars["importtime1"]+25],ages,vl_pars,virus2_inc+virus1_inc)

## Calculate distribution skew based on PDF
skew <- function(values,weights) {
  weights_std <- weights/sum(weights, na.rm=TRUE)
  xbar <- sum(values*weights_std, na.rm=TRUE)
  xi_xbar <- values - xbar
  return((sum(weights_std*xi_xbar^3))/((sum(weights_std*xi_xbar^2))^(3/2)))
}

## Summaries of infection ages over time:
age_res_std_overall <- age_res_overall/apply(age_res_overall, 1, sum, na.rm=TRUE)
age_res_std1 <- age_res1/apply(age_res1, 1, sum, na.rm=TRUE)
age_res_std2 <- age_res2/apply(age_res2, 1, sum, na.rm=TRUE)
age_mean1 <- tibble(t=times,mean_age=apply(age_res_std1, 1, function(res) sum(res*(1:lastday), na.rm=TRUE)),virus="Original variant")
age_mean2 <- tibble(t=times,mean_age=apply(age_res_std2, 1, function(res) sum(res*(1:lastday), na.rm=TRUE)),virus="New variant")
age_mean_overall <- tibble(t=times,mean_age=apply(age_res_std_overall, 1, function(res) sum(res*(1:lastday), na.rm=TRUE)),virus="Overall")
ages_deterministic_dat <- bind_rows(age_mean1,age_mean2,age_mean_overall)


## Summaries of strain 1 Ct values
dat1_skew <- dat1 %>% filter(ct < 40) %>% group_by(t) %>% summarize(skew1=skew(ct,density)) %>% mutate(virus="Original variant")
dat1_dist <- dat1 %>% filter(ct < 40) %>% group_by(t) %>%   mutate(density_scaled=density/sum(density)) %>%  mutate(cumu_density=cumsum(density_scaled)) 
dat1_mean <- dat1_dist %>% summarize(mean=sum(ct*density_scaled)) %>%  mutate(virus="Original variant")
dat1_med <- dat1_dist %>% filter(cumu_density >= 0.5) %>% filter(ct == min(ct)) %>%  mutate(virus="Original variant") %>% dplyr::select(ct, t, virus) %>% rename(median=ct)
dat1_low25 <- dat1_dist %>% filter(cumu_density >= 0.25) %>% filter(ct == min(ct)) %>%  mutate(virus="Original variant") %>% dplyr::select(ct, t, virus) %>% rename(low25=ct)
dat1_upp75 <- dat1_dist %>% filter(cumu_density >= 0.75) %>% filter(ct == min(ct)) %>%  mutate(virus="Original variant") %>% dplyr::select(ct, t, virus) %>% rename(upp75=ct)
dat1_summary <- left_join(dat1_med,dat1_low25) %>% left_join(dat1_upp75)

## Summaries of strain 2 Ct values
dat2_skew <- dat2 %>% filter(ct < 40) %>% group_by(t) %>% summarize(skew1=skew(ct,density)) %>% mutate(virus="New variant")
dat2_dist <- dat2 %>% filter(ct < 40) %>% group_by(t) %>%   mutate(density_scaled=density/sum(density)) %>%  mutate(cumu_density=cumsum(density_scaled)) 
dat2_mean <- dat2_dist %>% summarize(mean=sum(ct*density_scaled)) %>%  mutate(virus="New variant")
dat2_med <- dat2_dist %>% filter(cumu_density >= 0.5) %>% filter(ct == min(ct)) %>%  mutate(virus="New variant") %>% dplyr::select(ct, t, virus) %>% rename(median=ct)
dat2_low25 <- dat2_dist %>% filter(cumu_density >= 0.25) %>% filter(ct == min(ct)) %>%  mutate(virus="New variant") %>% dplyr::select(ct, t, virus) %>% rename(low25=ct)
dat2_upp75 <- dat2_dist %>% filter(cumu_density >= 0.75) %>% filter(ct == min(ct)) %>%  mutate(virus="New variant") %>% dplyr::select(ct, t, virus) %>% rename(upp75=ct)
dat2_summary <- left_join(dat2_med,dat2_low25) %>% left_join(dat2_upp75)

## Summaries of overall Ct values
dat_overall_skew <- dat_overall %>% filter(ct < 40) %>% group_by(t) %>% summarize(skew1=skew(ct,density)) %>% mutate(virus="Overall")
dat_overall_dist <- dat_overall %>% filter(ct < 40) %>% group_by(t) %>%   mutate(density_scaled=density/sum(density)) %>%  mutate(cumu_density=cumsum(density_scaled)) 
dat_overall_mean <- dat_overall_dist %>% summarize(mean=sum(ct*density_scaled)) %>%  mutate(virus="Overall")
dat_overall_med <- dat_overall_dist %>% filter(cumu_density >= 0.5) %>% filter(ct == min(ct)) %>%  mutate(virus="Overall") %>% dplyr::select(ct, t, virus) %>% rename(median=ct)
dat_overall_low25 <- dat_overall_dist %>% filter(cumu_density >= 0.25) %>% filter(ct == min(ct)) %>%  mutate(virus="Overall") %>% dplyr::select(ct, t, virus) %>% rename(low25=ct)
dat_overall_upp75 <- dat_overall_dist %>% filter(cumu_density >= 0.75) %>% filter(ct == min(ct)) %>%  mutate(virus="Overall") %>% dplyr::select(ct, t, virus) %>% rename(upp75=ct)
dat_overall_summary <- left_join(dat_overall_med,dat_overall_low25) %>% left_join(dat_overall_upp75)

## Combine summaries
combined_means <- bind_rows(dat1_mean,dat2_mean,dat_overall_mean)
combined_summaries <- bind_rows(dat1_summary,dat2_summary,dat_overall_summary)
combined_skews <- bind_rows(dat1_skew,dat2_skew,dat_overall_skew)
  
p_medians <- ggplot(combined_summaries) + 
  geom_ribbon(aes(x=t,ymin=low25,ymax=upp75,fill=virus),alpha=0.05) +
  geom_line(aes(x=t,y=median,col=virus,linetype=virus)) +
  geom_line(aes(x=t,y=low25,col=virus,linetype=virus)) +
  geom_line(aes(x=t,y=upp75,col=virus,linetype=virus)) +
  scale_fill_manual(values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_color_manual(values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_linetype_manual(values=c("Original variant"="solid","New variant"="solid","Overall"="dashed")) +
  scale_y_continuous(trans="reverse") +
  scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
  ylab("Median and 75% quantiles \nof detectable Ct values") +
  theme_overall +
  theme(legend.position="none") +theme_nice_axes + theme_no_x_axis +
  labs(tag="B")
  
p_skews <- ggplot(combined_skews %>% drop_na()) + 
  geom_line(aes(x=t,y=skew1,col=virus,linetype=virus))+
  scale_color_manual(values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_linetype_manual(values=c("Original variant"="solid","New variant"="solid","Overall"="dashed")) +
  scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
  #scale_y_continuous(expand=c(0,0)) +
  ylab("Skew of detectable Ct values") +
  theme_overall +
  xlab("Time") +
  theme(legend.position="none") + theme_nice_axes +
  labs(tag="C")

fig1 <- p_inc/p_medians/p_skews


###############################################
## 3) STOCHASTIC SIMULATION OF CT VALUES
###############################################
## Simulate a population of 1,000,000
population_n <- 1000000
## Simulate complete line list for individuals infection with the original or new variant
v1_linelist <- virosolver::simulate_observations_wrapper(virus1_inc*population_n,times=times,population_n=population_n) %>% 
  mutate(virus="Original variant") %>% filter(!is.na(infection_time)) %>% mutate(i=1:n())
v2_linelist <- virosolver::simulate_observations_wrapper(virus2_inc*population_n,times=times,population_n=population_n) %>% 
  mutate(virus="New variant") %>% filter(!is.na(infection_time))%>% mutate(i=1:n())
v2_linelist$i <- v2_linelist$i + max(v1_linelist$i)
complete_linelist <- bind_rows(v1_linelist, v2_linelist)
uninfected_linelist <- tibble(i=max(complete_linelist$i):population_n) %>%
                        mutate(infection_time=NA,is_infected=0,is_symp=0,incu_period=NA,onset_time=NA,
                              confirmation_delay=extraDistr::rdgamma(n(),5,2))
complete_linelist <- bind_rows(complete_linelist,uninfected_linelist)


reporting_prob <- 1
## Simulate situation where all individuals are observed at some point
observed_linelist <- simulate_reporting(complete_linelist %>% filter(is_infected == 1), 
                                        frac_report=reporting_prob,
                                        timevarying_prob=NULL,
                                        solve_times=times, 
                                        symptomatic=FALSE)



## Simulate Ct values from random cross sections
simulated_viral_loads <- simulate_viral_loads_wrapper(observed_linelist$sampled_individuals %>% filter(!is.na(infection_time)),kinetics_pars=vl_pars)

## Check incidence from linelist
complete_linelist %>% drop_na() %>% group_by(infection_time, virus) %>% summarize(inc=n()) %>% 
  ggplot() + geom_line(aes(x=infection_time,y=inc,col=as.factor(virus)))


p_inc2 <- ggplot(complete_linelist %>% drop_na() %>% group_by(infection_time, virus) %>% summarize(inc=n())) + 
  geom_line(aes(x=infection_time,y=inc,col=virus)) +
  geom_vline(xintercept=c(180),linetype="dashed",col="red") +
  geom_vline(xintercept=c(0),linetype="dashed",col="blue") +
  scale_color_manual(name="Variant",values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  xlab("Time") +
  scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
  #scale_y_continuous(expand=c(0,0)) +
  ylab("Incidence") +
  theme_overall + 
  theme(legend.position=c(0.8,0.8)) +
  theme_nice_axes+ theme_no_x_axis  +
  labs(tag="A")


var1_means <- simulated_viral_loads %>% filter(ct_obs < 40,virus=="Original variant") %>% group_by(sampled_time) %>% summarize(mean=mean(ct_obs),n=n())
var2_means <- simulated_viral_loads %>% filter(ct_obs < 40,virus=="New variant") %>% group_by(sampled_time) %>% summarize(mean=mean(ct_obs),n=n())
comb_means <- simulated_viral_loads %>% filter(ct_obs < 40) %>% group_by(sampled_time) %>% summarize(mean=mean(ct_obs),n=n())
n_min <- 10

p2 <- ggplot() + 
  geom_smooth(data=simulated_viral_loads %>% filter(ct_obs < 40, virus=="Original variant"),aes(x=sampled_time,y=ct_obs,col="Original variant",fill="Original variant"),alpha=0.1) +
  geom_line(data=var1_means %>% filter(n > n_min),aes(x=sampled_time,y=mean,col="Original variant"),size=0.1) +
  geom_smooth(data=simulated_viral_loads %>% filter(ct_obs < 40, virus=="New variant"),aes(x=sampled_time,y=ct_obs,col="New variant",fill="New variant"),alpha=0.1) +
  geom_line(data=var2_means %>% filter(n > n_min),aes(x=sampled_time,y=mean,col="New variant"),size=0.1) +
  geom_smooth(data=simulated_viral_loads %>% filter(ct_obs < 40),aes(x=sampled_time,y=ct_obs,col="Overall",fill="Overall"),alpha=0.1,linetype="dashed") +
  geom_line(data=comb_means %>% filter(n > n_min),aes(x=sampled_time,y=mean,col="Overall"),size=0.1) +
  scale_color_manual(values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_fill_manual(values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_y_continuous(trans="reverse") +
  scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
  theme_overall +
  theme(legend.position="none") +theme_nice_axes + theme_no_x_axis +
  xlab("Time") +
  ylab("Smoothed detectable Ct values") +
  labs(tag="B")
  

var1_skews <- simulated_viral_loads %>% filter(ct_obs < 40,virus=="New variant") %>% group_by(sampled_time) %>% summarize(skewness=moments::skewness(ct_obs),n=n())
var2_skews <- simulated_viral_loads %>% filter(ct_obs < 40,virus=="Original variant") %>% group_by(sampled_time) %>% summarize(skewness=moments::skewness(ct_obs),n=n())
comb_skews <- simulated_viral_loads %>% filter(ct_obs < 40) %>% group_by(sampled_time) %>% summarize(skewness=moments::skewness(ct_obs),n=n())

p3 <- ggplot()+
  geom_smooth(data=var1_skews %>% filter(n > n_min),aes(x=sampled_time,y=skewness,col="New variant",fill="New variant"),alpha=0.25) +
  geom_line(data=var1_skews %>% filter(n > n_min),aes(x=sampled_time,y=skewness,col="New variant"),size=0.1) +
  geom_smooth(data=var2_skews %>% filter(n > n_min), aes(x=sampled_time,y=skewness,col="Original variant",fill="Original variant"),alpha=0.25)+
  geom_line(data=var2_skews %>% filter(n > n_min), aes(x=sampled_time,y=skewness,col="Original variant"),size=0.1)+
  geom_smooth(data=comb_skews %>% filter(n > n_min),aes(x=sampled_time,y=skewness,col="Overall",fill="Overall"),alpha=0.25,linetype="dashed") +
  geom_line(data=comb_skews %>% filter(n > n_min), aes(x=sampled_time,y=skewness,col="Overall"),size=0.1) +
  theme(legend.position="none",legend.title=element_blank()) +
  scale_color_manual(values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_fill_manual(values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
  theme_overall +
  theme(legend.position="none") +theme_nice_axes + 
  xlab("Time") +
  ylab("Skewness of detectable Ct values")+
  labs(tag="C")
  
fig2 <- p_inc2/p2/p3
dat_inf_times_by_variant <- simulated_viral_loads %>% filter(ct_obs < 40) %>% 
  group_by(sampled_time,virus) %>% 
  summarize(mean_inf_time=mean(days_since_infection),n=n()) %>%
  filter(n > n_min)
dat_inf_times_combined <- simulated_viral_loads %>% filter(ct_obs < 40) %>% 
  group_by(sampled_time) %>% 
  summarize(mean_inf_time=mean(days_since_infection),n=n()) %>%
  filter(n > n_min) %>% mutate(virus="Overall")
dat_inf_times <- bind_rows(dat_inf_times_by_variant, dat_inf_times_combined)

fig3a <- ggplot(ages_deterministic_dat) +
  geom_line(aes(x=t,y=mean_age,col=virus),size=1) +
  scale_color_manual(name="Variant",values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
  scale_y_continuous(breaks=seq(0,25,by=5)) +
  coord_cartesian(ylim=c(5,25)) +
  theme_overall +
  theme(legend.position=c(0.25,0.75)) +theme_nice_axes + 
  xlab("Time") +
  ylab("Mean time since infection of sampled individuals") +
  labs(tag="A")
fig3b <- ggplot(dat_inf_times) +
  geom_line(aes(x=sampled_time,y=mean_inf_time,col=virus),size=0.1) +
  geom_smooth(data=simulated_viral_loads %>% filter(ct_obs < 40), aes(x=sampled_time,y=days_since_infection,col=virus,fill=virus),alpha=0.25) +
  geom_smooth(data=simulated_viral_loads %>% filter(ct_obs < 40), aes(x=sampled_time,y=days_since_infection,col="Overall",fill="Overall"),alpha=0.25) +
  scale_color_manual(name="Variant",values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_fill_manual(name="Variant",values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
  scale_y_continuous(breaks=seq(0,25,by=5)) +
  coord_cartesian(ylim=c(5,25)) +
  theme_overall +
  theme(legend.position=c(0.25,0.75)) +theme_nice_axes + 
  xlab("Time") +
  ylab("Mean time since infection of sampled individuals") +
  labs(tag="B")

fig3 <- fig3a/fig3b

###############################################
## 4) STOCHASTIC SIMULATION OF CT VALUES, SYMPTOMATIC REPORTING
###############################################
reporting_prob <- 1
## Simulate situation where all individuals are observed at some point
observed_linelist_symptomatic <- simulate_reporting(complete_linelist %>% filter(is_infected == 1), 
                                        frac_report=reporting_prob,
                                        timevarying_prob=NULL,
                                        solve_times=times, 
                                        symptomatic=TRUE)


## Simulate Ct values from random cross sections
simulated_viral_loads_symptomatic <- simulate_viral_loads_wrapper(observed_linelist_symptomatic$sampled_individuals %>% filter(!is.na(infection_time)),kinetics_pars=vl_pars)

var1_means <- simulated_viral_loads_symptomatic %>% filter(ct_obs < 40,virus=="Original variant") %>% group_by(sampled_time) %>% summarize(mean=mean(ct_obs),n=n())
var2_means <- simulated_viral_loads_symptomatic %>% filter(ct_obs < 40,virus=="New variant") %>% group_by(sampled_time) %>% summarize(mean=mean(ct_obs),n=n())
comb_means <- simulated_viral_loads_symptomatic %>% filter(ct_obs < 40) %>% group_by(sampled_time) %>% summarize(mean=mean(ct_obs),n=n())
n_min <- 10

p2_sympt <- ggplot() + 
  geom_smooth(data=simulated_viral_loads_symptomatic %>% filter(ct_obs < 40, virus=="Original variant"),aes(x=sampled_time,y=ct_obs,col="Original variant",fill="Original variant"),alpha=0.1) +
  #geom_line(data=var1_means %>% filter(n > n_min),aes(x=sampled_time,y=mean,col="Original variant"),size=0.1) +
  geom_smooth(data=simulated_viral_loads_symptomatic %>% filter(ct_obs < 40, virus=="New variant"),aes(x=sampled_time,y=ct_obs,col="New variant",fill="New variant"),alpha=0.1) +
  #geom_line(data=var2_means %>% filter(n > n_min),aes(x=sampled_time,y=mean,col="New variant"),size=0.1) +
  geom_smooth(data=simulated_viral_loads_symptomatic %>% filter(ct_obs < 40),aes(x=sampled_time,y=ct_obs,col="Overall",fill="Overall"),alpha=0.1,linetype="dashed") +
  #geom_line(data=comb_means %>% filter(n > n_min),aes(x=sampled_time,y=mean,col="Overall"),size=0.1) +
  scale_color_manual(values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_fill_manual(values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_y_continuous(trans="reverse") +
  scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
  theme_overall +
  theme(legend.position="none") +theme_nice_axes + theme_no_x_axis +
  xlab("Time") +
  ylab("Smoothed detectable Ct values") +
  labs(tag="B")


var1_skews <- simulated_viral_loads_symptomatic %>% filter(ct_obs < 40,virus=="New variant") %>% group_by(sampled_time) %>% summarize(skewness=moments::skewness(ct_obs),n=n())
var2_skews <- simulated_viral_loads_symptomatic %>% filter(ct_obs < 40,virus=="Original variant") %>% group_by(sampled_time) %>% summarize(skewness=moments::skewness(ct_obs),n=n())
comb_skews <- simulated_viral_loads_symptomatic %>% filter(ct_obs < 40) %>% group_by(sampled_time) %>% summarize(skewness=moments::skewness(ct_obs),n=n())

p3_sympt <- ggplot()+
  geom_smooth(data=var1_skews %>% filter(n > n_min),aes(x=sampled_time,y=skewness,col="New variant",fill="New variant"),alpha=0.25) +
  #geom_line(data=var1_skews %>% filter(n > n_min),aes(x=sampled_time,y=skewness,col="New variant"),size=0.1) +
  geom_smooth(data=var2_skews %>% filter(n > n_min), aes(x=sampled_time,y=skewness,col="Original variant",fill="Original variant"),alpha=0.25)+
  #geom_line(data=var2_skews %>% filter(n > n_min), aes(x=sampled_time,y=skewness,col="Original variant"),size=0.1)+
  geom_smooth(data=comb_skews %>% filter(n > n_min),aes(x=sampled_time,y=skewness,col="Overall",fill="Overall"),alpha=0.25,linetype="dashed") +
  #geom_line(data=comb_skews %>% filter(n > n_min), aes(x=sampled_time,y=skewness,col="Overall"),size=0.1) +
  theme(legend.position="none",legend.title=element_blank()) +
  scale_color_manual(values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_fill_manual(values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
  theme_overall +
  theme(legend.position="none") +theme_nice_axes + 
  xlab("Time") +
  ylab("Skewness of detectable Ct values")+
  labs(tag="C")

fig4 <- p_inc2/p2_sympt/p3_sympt

###############################################
## 6) STOCHASTIC SIMULATION OF CT VALUES, DIFFERENT KINETICS
###############################################
vl_pars_var1 <- vl_pars_var2 <- vl_pars
vl_pars_var2["viral_peak"] <- 17
vl_pars_var2["t_switch"] <- 15

## Simulate Ct values from random cross sections
simulated_viral_loads_var1 <- simulate_viral_loads_wrapper(observed_linelist$sampled_individuals %>% filter(!is.na(infection_time),virus=="Original variant"),kinetics_pars=vl_pars_var1)
simulated_viral_loads_var2 <- simulate_viral_loads_wrapper(observed_linelist$sampled_individuals %>% filter(!is.na(infection_time),virus=="New variant"),kinetics_pars=vl_pars_var2)
simulated_viral_loads <- bind_rows(simulated_viral_loads_var1,simulated_viral_loads_var2)

var1_means <- simulated_viral_loads %>% filter(ct_obs < 40,virus=="Original variant") %>% group_by(sampled_time) %>% summarize(mean=mean(ct_obs),n=n())
var2_means <- simulated_viral_loads %>% filter(ct_obs < 40,virus=="New variant") %>% group_by(sampled_time) %>% summarize(mean=mean(ct_obs),n=n())
comb_means <- simulated_viral_loads %>% filter(ct_obs < 40) %>% group_by(sampled_time) %>% summarize(mean=mean(ct_obs),n=n())
n_min <- 10

p2_6 <- ggplot() + 
  geom_smooth(data=simulated_viral_loads %>% filter(ct_obs < 40, virus=="Original variant"),aes(x=sampled_time,y=ct_obs,col="Original variant",fill="Original variant"),alpha=0.1) +
  geom_line(data=var1_means %>% filter(n > n_min),aes(x=sampled_time,y=mean,col="Original variant"),size=0.1) +
  geom_smooth(data=simulated_viral_loads %>% filter(ct_obs < 40, virus=="New variant"),aes(x=sampled_time,y=ct_obs,col="New variant",fill="New variant"),alpha=0.1) +
  geom_line(data=var2_means %>% filter(n > n_min),aes(x=sampled_time,y=mean,col="New variant"),size=0.1) +
  geom_smooth(data=simulated_viral_loads %>% filter(ct_obs < 40),aes(x=sampled_time,y=ct_obs,col="Overall",fill="Overall"),alpha=0.1,linetype="dashed") +
  geom_line(data=comb_means %>% filter(n > n_min),aes(x=sampled_time,y=mean,col="Overall"),size=0.1) +
  scale_color_manual(values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_fill_manual(values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_y_continuous(trans="reverse") +
  scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
  theme_overall +
  theme(legend.position="none") +theme_nice_axes + theme_no_x_axis +
  xlab("Time") +
  ylab("Smoothed detectable Ct values") +
  labs(tag="B")


var1_skews <- simulated_viral_loads %>% filter(ct_obs < 40,virus=="New variant") %>% group_by(sampled_time) %>% summarize(skewness=moments::skewness(ct_obs),n=n())
var2_skews <- simulated_viral_loads %>% filter(ct_obs < 40,virus=="Original variant") %>% group_by(sampled_time) %>% summarize(skewness=moments::skewness(ct_obs),n=n())
comb_skews <- simulated_viral_loads %>% filter(ct_obs < 40) %>% group_by(sampled_time) %>% summarize(skewness=moments::skewness(ct_obs),n=n())

p3_6 <- ggplot()+
  geom_smooth(data=var1_skews %>% filter(n > n_min),aes(x=sampled_time,y=skewness,col="New variant",fill="New variant"),alpha=0.25) +
  geom_line(data=var1_skews %>% filter(n > n_min),aes(x=sampled_time,y=skewness,col="New variant"),size=0.1) +
  geom_smooth(data=var2_skews %>% filter(n > n_min), aes(x=sampled_time,y=skewness,col="Original variant",fill="Original variant"),alpha=0.25)+
  geom_line(data=var2_skews %>% filter(n > n_min), aes(x=sampled_time,y=skewness,col="Original variant"),size=0.1)+
  geom_smooth(data=comb_skews %>% filter(n > n_min),aes(x=sampled_time,y=skewness,col="Overall",fill="Overall"),alpha=0.25,linetype="dashed") +
  geom_line(data=comb_skews %>% filter(n > n_min), aes(x=sampled_time,y=skewness,col="Overall"),size=0.1) +
  theme(legend.position="none",legend.title=element_blank()) +
  scale_color_manual(values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_fill_manual(values=c("Original variant"="blue","New variant"="red","Overall"="black")) +
  scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
  theme_overall +
  theme(legend.position="none") +theme_nice_axes + 
  xlab("Time") +
  ylab("Skewness of detectable Ct values")+
  labs(tag="C")

fig7 <- p_inc2/p2_6/p3_6

###############################################
## SAVE PLOTS
###############################################
ggsave("figures/fig1.png",fig1,height=7,width=7,dpi=300,units = "in")
ggsave("figures/fig2.png",fig2,height=7,width=7,dpi=300,units = "in")
ggsave("figures/fig3.png",fig3,height=7,width=7,dpi=300,units = "in")
ggsave("figures/fig4.png",fig4,height=7,width=7,dpi=300,units = "in")
ggsave("figures/fig7.png",fig7,height=7,width=7,dpi=300,units = "in")

