library(tidyverse)
library(jahR)

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

samp_time <- 50
tmp1 <- (virus1_inc[seq(samp_time-35,samp_time,by=1)])
early_mean_gr_virus1 <- mean(log(tmp1[2:length(tmp1)]/tmp1[1:(length(tmp1)-1)]))
tmp2 <- (virus2_inc[seq(samp_time-35,samp_time,by=1)])
early_mean_gr_virus2 <- mean(log(tmp2[2:length(tmp2)]/tmp2[1:(length(tmp2)-1)]))

early_grs <- tibble(name=c("Original variant","New variant","Difference"),
                    gr=c(early_mean_gr_virus1, early_mean_gr_virus2, early_mean_gr_virus2-early_mean_gr_virus1))
early_grs$name <- factor(early_grs$name,levels=c("Original variant","New variant","Difference"))

pars["importtime2"] <- 180
seir_dynamics <- run_2strain_seir_simulation(pars, states,times)
virus1_inc <- seir_dynamics$virus1_inc
virus2_inc <- seir_dynamics$virus2_inc
virus_inc <- virus1_inc + virus2_inc

samp_time <- 270
tmp1 <- (virus1_inc[seq(samp_time-35,samp_time,by=1)])
late_mean_gr_virus1 <- mean(log(tmp1[2:length(tmp1)]/tmp1[1:(length(tmp1)-1)]))
tmp2 <- (virus2_inc[seq(samp_time-35,samp_time,by=1)])
late_mean_gr_virus2 <- mean(log(tmp2[2:length(tmp2)]/tmp2[1:(length(tmp2)-1)]))
late_grs <- tibble(name=c("Original variant","New variant","Difference"),
                    gr=c(late_mean_gr_virus1, late_mean_gr_virus2, late_mean_gr_virus2-late_mean_gr_virus1))
late_grs$name <- factor(late_grs$name,levels=c("Original variant","New variant","Difference"))


## Create plot of growth rate estimates for day 270
setwd("~/Downloads/results/virosolver_late/growth_rates")
files <- list.files()
all_res <- NULL

for(i in seq_along(files)){
  all_res[[i]] <- read_csv(files[i])
}
all_res <- do.call("bind_rows", all_res)

all_res$name <- factor(all_res$name, levels=c("Original variant","New variant", "Difference"))
all_res$samp_size <- paste0(all_res$samp_size," samples")
all_res$samp_size <- factor(all_res$samp_size, levels=c("25 samples","50 samples","100 samples","250 samples","500 samples"))

summaries <- all_res %>% 
  group_by(samp_size,run_ID,name) %>% 
  summarize(mean_est=mean(value),
           lower_95=quantile(value,0.025),
           lower_50=quantile(value,0.25),
           median_est=median(value),
           upper_50=quantile(value,0.75),
           upper_95=quantile(value,0.975))

p1 <- ggplot(summaries) + geom_jitter(aes(x=name,y=mean_est),height=0,width=0.25) + facet_wrap(~samp_size) +
  scale_y_continuous(limits=c(-0.2,0.2)) +
  geom_hline(yintercept=0,col="red",linetype="dashed")
fig_late <- ggplot(summaries %>% group_by(samp_size,name) %>% mutate(id=1:n())) +
  geom_hline(yintercept=0,size=0.5,col="grey40") +
  geom_errorbar(aes(x=id,ymin=lower_95,ymax=upper_95,col=name),size=0.3,width=0,alpha=0.5) +
  geom_errorbar(aes(x=id,ymin=lower_50,ymax=upper_50,col=name),size=0.3,width=0) +
  geom_point(aes(x=id,y=mean_est,col=name),size=0.75,shape=15) + 
  geom_hline(data=late_grs,aes(yintercept=gr),col="red",linetype="dashed") +
  variant_color_scale +
  facet_grid(samp_size~name)+
  ylab("Estimated growth rate") +
  xlab("Simulation ID") +
  scale_y_continuous(limits=c(-0.25,0.25)) +
  theme_overall()+
  theme(legend.position="none", 
        strip.background = element_blank(),
        strip.text=element_text(size=10,face="bold")) 


## Create plot of growth rate estimates for day 50
setwd("~/Downloads/results/virosolver_diff_directions/growth_rates")
files <- list.files()
all_res <- NULL

for(i in seq_along(files)){
  all_res[[i]] <- read_csv(files[i])
}
all_res <- do.call("bind_rows", all_res)

all_res$name <- factor(all_res$name, levels=c("Original variant","New variant", "Difference"))
all_res$samp_size <- paste0(all_res$samp_size," samples")
all_res$samp_size <- factor(all_res$samp_size, levels=c("25 samples","50 samples","100 samples","250 samples","500 samples"))

summaries <- all_res %>% 
  group_by(samp_size,run_ID,name) %>% 
  summarize(mean_est=mean(value),
            lower_95=quantile(value,0.025),
            lower_50=quantile(value,0.25),
            median_est=median(value),
            upper_50=quantile(value,0.75),
            upper_95=quantile(value,0.975))

p1 <- ggplot(summaries) + geom_jitter(aes(x=name,y=mean_est),height=0,width=0.25) + facet_wrap(~samp_size) +
  scale_y_continuous(limits=c(-0.2,0.2)) +
  geom_hline(yintercept=0,col="red",linetype="dashed")
fig_early <- ggplot(summaries %>% group_by(samp_size,name) %>% mutate(id=1:n())) +
  geom_hline(yintercept=0,size=0.5,col="grey40") +
  geom_errorbar(aes(x=id,ymin=lower_95,ymax=upper_95,col=name),size=0.3,width=0,alpha=0.5) +
  geom_errorbar(aes(x=id,ymin=lower_50,ymax=upper_50,col=name),size=0.3,width=0) +
  geom_point(aes(x=id,y=mean_est,col=name),size=0.75,shape=15) + 
  geom_hline(data=early_grs,aes(yintercept=gr),col="red",linetype="dashed") +
  variant_color_scale +
  facet_grid(samp_size~name)+
  ylab("Estimated growth rate") +
  xlab("Simulation ID") +
  scale_y_continuous(limits=c(-0.25,0.25)) +
  theme_overall()+
  theme(legend.position="none", 
        strip.background = element_blank(),
        strip.text=element_text(size=10,face="bold")) 



## Plot viral kinetics difference estimates

vl_pars1 <- vl_pars
vl_pars2 <- vl_pars_both
real_scales <- tibble(name=c("Relative peak\n Ct value","Relative duration\n of initial waning"),value=c(vl_pars2["viral_peak"]/vl_pars1["viral_peak"], vl_pars2["t_switch"]/vl_pars1["t_switch"]))


setwd("~/Downloads/results/virosolver_diff_directions/viral_kinetics/")
files <- list.files()
all_res <- NULL
for(i in seq_along(files)){
  all_res[[i]] <- read_csv(files[i])
}
all_res <- do.call("bind_rows", all_res)
all_res$samp_size <- paste0(all_res$samp_size," samples")
all_res$samp_size <- factor(all_res$samp_size, levels=c("25 samples","50 samples","100 samples","250 samples","500 samples"))

summaries <- all_res %>% 
  group_by(samp_size,run_ID,name) %>% 
  summarize(mean_est=mean(value),
            lower_95=quantile(value,0.025),
            lower_50=quantile(value,0.25),
            median_est=median(value),
            upper_50=quantile(value,0.75),
            upper_95=quantile(value,0.975))
fig_early_viral <- ggplot(summaries %>% group_by(samp_size) %>% mutate(id=1:n())) +
  geom_hline(yintercept=1,size=0.5,col="grey40") +
  geom_errorbar(aes(x=id,ymin=lower_95,ymax=upper_95,col=name),size=0.3,width=0,alpha=0.5) +
  geom_errorbar(aes(x=id,ymin=lower_50,ymax=upper_50,col=name),size=0.3,width=0) +
  geom_point(aes(x=id,y=mean_est,col=name),size=0.75,shape=15) + 
  geom_hline(data=real_scales,aes(yintercept=value),linetype="dashed",col="red") +
  variant_color_scale +
  facet_grid(samp_size~name)+
  ylab("Estimate") +
  xlab("Simulation ID") +
  theme_overall()+
  theme(legend.position="none", 
        strip.background = element_blank(),
        strip.text=element_text(size=10,face="bold")) 


## Plot viral kinetics difference estimates

vl_pars1 <- vl_pars
vl_pars2 <- vl_pars_both
real_scales <- tibble(name=c("Relative peak\n Ct value","Relative duration\n of initial waning"),value=c(vl_pars2["viral_peak"]/vl_pars1["viral_peak"], vl_pars2["t_switch"]/vl_pars1["t_switch"]))


setwd("~/Downloads/results/virosolver_late/viral_kinetics/")
files <- list.files()
all_res <- NULL
for(i in seq_along(files)){
  all_res[[i]] <- read_csv(files[i])
}
all_res <- do.call("bind_rows", all_res)
all_res$samp_size <- paste0(all_res$samp_size," samples")
all_res$samp_size <- factor(all_res$samp_size, levels=c("25 samples","50 samples","100 samples","250 samples","500 samples"))

summaries <- all_res %>% 
  group_by(samp_size,run_ID,name) %>% 
  summarize(mean_est=mean(value),
            lower_95=quantile(value,0.025),
            lower_50=quantile(value,0.25),
            median_est=median(value),
            upper_50=quantile(value,0.75),
            upper_95=quantile(value,0.975))
fig_late_viral <- ggplot(summaries %>% group_by(samp_size) %>% mutate(id=1:n())) +
  geom_hline(yintercept=1,size=0.5,col="grey40") +
  geom_errorbar(aes(x=id,ymin=lower_95,ymax=upper_95,col=name),size=0.3,width=0,alpha=0.5) +
  geom_errorbar(aes(x=id,ymin=lower_50,ymax=upper_50,col=name),size=0.3,width=0) +
  geom_point(aes(x=id,y=mean_est,col=name),size=0.75,shape=15) + 
  geom_hline(data=real_scales,aes(yintercept=value),linetype="dashed",col="red") +
  variant_color_scale +
  facet_grid(samp_size~name)+
  ylab("Estimate") +
  xlab("Simulation ID") +
  theme_overall()+
  theme(legend.position="none", 
        strip.background = element_blank(),
        strip.text=element_text(size=10,face="bold")) 


setwd("~/Documents/GitHub/variant_viral_loads/")
ggsave("figures/virosolver_compare_early.pdf",fig_early,width=8,height=8)
ggsave("figures/virosolver_compare_early.png",fig_early,width=8,height=8,units="in",dpi=300)

ggsave("figures/virosolver_compare_late.pdf",fig_late,width=8,height=8)
ggsave("figures/virosolver_compare_late.png",fig_late,width=8,height=8,units="in",dpi=300)

ggsave("figures/virosolver_compare_early_viral.pdf",fig_early_viral,width=8,height=8)
ggsave("figures/virosolver_compare_early_viral.png",fig_early_viral,width=8,height=8,units="in",dpi=300)

ggsave("figures/virosolver_compare_late_viral.pdf",fig_late_viral,width=8,height=8)
ggsave("figures/virosolver_compare_late_viral.png",fig_late_viral,width=8,height=8,units="in",dpi=300)




