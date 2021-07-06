###############################################
## DISTINGUISHING GROWTH RATES AND VIRAL KINETICS FROM VARIANTS
## James Hay
## jhay@hsph.harvard.edu
## June 8, 2021
###############################################
library(virosolver)
library(tidyverse)
library(patchwork)
library(extraDistr)
library(virosolver)
library(ggpubr)
#library(lazymcmc)
devtools::load_all("~/Documents/GitHub/lazymcmc/")
library(rethinking)
## Where to perform the simulations
HOME_WD <- "~"
HOME_WD <- "~/Documents/GitHub/variant_viral_loads/"
setwd(HOME_WD)


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
          importtime2 = 0, ## Time of importation of strain 2
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

samp_time <- 50
vl_pars1 <- vl_pars
vl_pars2 <- vl_pars_both
cts_1 <- tibble(ct=simulate_cross_section(vl_pars1, ages, virus1_inc,obs_time=samp_time,N=100),variant="Original variant")
cts_2 <- tibble(ct=simulate_cross_section(vl_pars2, ages, virus2_inc,obs_time=samp_time,N=100),variant="New variant, different kinetics")
cts_sim_comb <- bind_rows(cts_1,cts_2) %>% mutate(t=35)
cts_sim_comb$virus <- factor(cts_sim_comb$variant,levels=variant_levels)
p_ct_samp <- ggplot(cts_sim_comb) + 
  geom_violin(aes(x=virus,y=ct,fill=virus),alpha=0.1,trim=FALSE,col=NA,width=0.75) +
  geom_dotplot(aes(x=virus,y=ct,fill=virus),binaxis="y",stackdir="center",binwidth=1,dotsize=1) + 
  geom_errorbar(data=cts_sim_comb%>%group_by(virus)%>%summarize(median_ct=median(ct)),
                aes(y=median_ct,ymin=median_ct,ymax=median_ct,x=virus),size=0.5,col="grey10") +
  variant_fill_scale + variant_color_scale +
  scale_y_continuous(trans="reverse",limits=c(40,7)) +
  ylab("Ct value") +
  theme_overall + theme_nice_axes + theme(legend.position="none",axis.title.x=element_blank(),plot.title=element_text(size=7)) +
  scale_x_discrete(labels = c("Original variant" = "Original variant",
                              "New variant, different kinetics"="New variant,\ndifferent kinetics")) 


tmp1 <- (virus1_inc[seq(samp_time-35,samp_time,by=1)])
mean_gr_virus1 <- mean(log(tmp1[2:length(tmp1)]/tmp1[1:(length(tmp1)-1)]))
tmp2 <- (virus2_inc[seq(samp_time-35,samp_time,by=1)])
mean_gr_virus2 <- mean(log(tmp2[2:length(tmp2)]/tmp2[1:(length(tmp2)-1)]))


grs %>% filter(t <= samp_time & t >= (samp_time-35)) %>% ggplot() + geom_line(aes(x=t,y=gr,col=virus))


virosolver_pars <- read.csv("pars/partab_exp_model_compare.csv")

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

#mcmc_pars <- c("iterations"=100000,"popt"=0.44,"opt_freq"=2000,
#               "thin"=10,"adaptive_period"=50000,"save_block"=1000)


n_temperatures <- 10
mcmcPars_ct_pt <- list("iterations"=50000,"popt"=0.234,"opt_freq"=1000,
                       "thin"=1,"adaptive_period"=50000,"save_block"=1000,
                       "temperature" = seq(1,101,length.out=n_temperatures),
                       "parallel_tempering_iter" = 5,"max_adaptive_period" = 50000, 
                       "adaptiveLeeway" = 0.2, "max_total_iterations" = 50000)

#virosolver_pars[!(virosolver_pars$names %in% c("beta","beta_alt","viral_peak_scale","t_switch_scale")),"fixed"] <- 1


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
                   filename="chains/example",
                   CREATE_POSTERIOR_FUNC=create_posterior_func_compare,
                   mvrPars=mvrPars,
                   OPT_TUNING=0.2,
                   use_pos=TRUE,
                   solve_likelihood=TRUE)

chain <- read.csv(output$file)
chain <- chain[chain$sampno > mcmcPars_ct_pt["adaptive_period"],]
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


chain_scales <- chain %>% dplyr::select(sampno, viral_peak_scale, t_switch_scale) %>% pivot_longer(-sampno)
scale_key <- c("viral_peak_scale"="Relative peak\n Ct value","t_switch_scale"="Relative duration\n of initial waning")
chain_scales$name <- scale_key[chain_scales$name]
real_scales <- tibble(name=scale_key,value=c(vl_pars2["viral_peak"]/vl_pars1["viral_peak"], vl_pars2["t_switch"]/vl_pars1["t_switch"]))
p2 <- ggplot(chain_scales) + 
  geom_hline(yintercept=1,linetype="dashed")+
  geom_violin(aes(x=name,y=value),alpha=0.5,draw_quantiles=c(0.025,0.5,0.975),fill="red") + 
  geom_point(data=real_scales,aes(x=name,y=value)) +
  #scale_y_continuous(limits=c(0,7)) +
  scale_y_log10() +
  theme_overall +
  ylab("Relative value of new variant") +
  xlab("")


real_v1_gr <- tibble(t=0:35,prob_infection=tmp1/sum(tmp1))
v1_preds <- virosolver::plot_prob_infection(chain, 1000,exponential_growth_model, 0:35,true_prob_infection = real_v1_gr)
v1_preds_quants <- v1_preds$predictions %>% group_by(t) %>% summarize(lower95=quantile(prob_infection,0.025),lower50=quantile(prob_infection,0.25),
                                                   median=median(prob_infection),
                                                   upper50=quantile(prob_infection,0.75),upper95=quantile(prob_infection,0.975)) 
p_v1_preds <- ggplot(v1_preds_quants) + 
  geom_ribbon(aes(x=t,ymin=lower95,ymax=upper95),fill="blue",alpha=0.1) +
  geom_ribbon(aes(x=t,ymin=lower50,ymax=upper50),fill="blue",alpha=0.5) +
  geom_line(aes(x=t,y=median),col="blue") +
  geom_line(data=real_v1_gr,aes(x=t,y=prob_infection),linetype="dashed",col="black",size=0.75) +
  theme_overall +
  ylab("Original variant growth rate") +
  theme_no_x_axis
  

chain1 <- chain
chain1$beta <- chain1$beta_alt
real_v2_gr <- tibble(t=0:35,prob_infection=tmp2/sum(tmp2))
v2_preds <- virosolver::plot_prob_infection(chain1, 1000,exponential_growth_model, 0:35,true_prob_infection = real_v2_gr)
v2_preds_quants <- v2_preds$predictions %>% group_by(t) %>% summarize(lower95=quantile(prob_infection,0.025),lower50=quantile(prob_infection,0.25),
                                                                      median=median(prob_infection),
                                                                      upper50=quantile(prob_infection,0.75),upper95=quantile(prob_infection,0.975)) 
p_v2_preds <- ggplot(v2_preds_quants) + 
  geom_ribbon(aes(x=t,ymin=lower95,ymax=upper95),fill="red",alpha=0.1) +
  geom_ribbon(aes(x=t,ymin=lower50,ymax=upper50),fill="red",alpha=0.5) +
  geom_line(aes(x=t,y=median),col="red") +
  geom_line(data=real_v2_gr,aes(x=t,y=prob_infection),linetype="dashed",col="black",size=0.75) +
  theme_overall+
  ylab("New variant growth rate") +
  xlab("Time (35 days prior to sample)")

v1_preds_quants$virus <- "Original variant"
v2_preds_quants$virus <- "New variant, different kinetics"

v_preds_comb <- bind_rows(v1_preds_quants,v2_preds_quants)
p_comb_preds <- ggplot(v_preds_comb) + 
  geom_ribbon(aes(x=t,ymin=lower95,ymax=upper95,fill=virus),alpha=0.1) +
  geom_ribbon(aes(x=t,ymin=lower50,ymax=upper50,fill=virus,col=virus),alpha=0.5,size=0.2,linetype="dashed") +
  geom_line(aes(x=t,y=median,col=virus)) +
  variant_color_scale + variant_fill_scale +
  #geom_line(data=real_v2_gr,aes(x=t,y=prob_infection),linetype="dashed",col="black",size=0.75) +
  theme_overall+
  theme_nice_axes +
  theme(legend.position="none") +
  scale_x_continuous(breaks=seq(0,35,by=5),labels=rev(0-seq(0,35,by=5))) +
  ylab("Relative probability of infection") +
  xlab("Days (relative to sample date)")

p_compare_kinetics <- p_compare_estimated_curves(chain,ages=0:35,N=100,nsamp=1000)


p_main2 <- (p_ct_samp+labs(tag="D"))/
              (p_compare_kinetics+labs(tag="E"))/
              (p_comb_preds+labs(tag="F"))

#p_main <- p1+p2+p_v1_preds+p_v2_preds + plot_layout(ncol=2,byrow=TRUE)
#p_main
#ggsave("figures/estimate_diffs.png",p_main,height=6,width=8,units="in",dpi=300)
