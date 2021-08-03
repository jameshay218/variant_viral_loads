theme_no_x_axis <- theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank(),axis.line.x=element_blank())
theme_nice_axes <- theme(axis.line.x=element_line(size=0.5,color="black"),axis.line.y=element_line(size=0.5,color="black"),panel.border = element_blank())
theme_overall <- theme_bw() + theme(axis.text=element_text(size=7),axis.title = element_text(size=8),plot.tag = element_text(size=10,face="bold"),
                                    legend.text=element_text(size=7),legend.title=element_text(size=8),strip.background = element_blank(),
                                    strip.text=element_text(size=8,face="bold"))

variant_levels <- c("Overall","Original variant","New variant","New variant, same kinetics","New variant, different kinetics")

color_key <- c("Overall"="black",
               "Original variant"="#0072B2",
               "New variant"="#D55E00",
               "New variant, same kinetics"="#009E73",
               "New variant, different kinetics"="#D55E00",
               "Difference"="gray40")

linetype_key <- c("Original variant"="solid",
                  "New variant, different kinetics"="solid",
                  "New variant"="solid",
                  "New variant, same kinetics"="solid",
                  "Overall"="dashed")
  
size_key <- c("Original variant"=0.5,
              "New variant, different kinetics"=0.75,
              "New variant"=0.5,
              "New variant, same kinetics"=0.5,
              "Overall"=0.5)

variant_color_scale <- scale_color_manual(name="Variant",values=color_key,drop=TRUE) 
variant_fill_scale <- scale_fill_manual(name="Variant",values=color_key,drop=TRUE) 

variant_color_scale_fig1 <- scale_color_manual(name="Variant",values=color_key[c("Original variant","New variant","Overall")],drop=TRUE) 
variant_color_scale_fig2 <- scale_color_manual(name="Variant",drop=TRUE,values=color_key[c("Original variant","New variant, same kinetics","New variant, different kinetics")]) 

variant_color_scale_min <- scale_color_manual(name="Variant",values=color_key[c("Original variant","New variant")],drop=TRUE) 
variant_fill_scale_min <- scale_fill_manual(name="Variant",values=color_key[c("Original variant","New variant")],drop=TRUE) 
variant_fill_scale_fig2 <- scale_fill_manual(name="Variant",drop=TRUE,values=color_key[c("Original variant","New variant, same kinetics","New variant, different kinetics")]) 

variant_linetype_scale <- scale_linetype_manual(name="Variant",values=linetype_key,drop=TRUE)
variant_size_scale <- scale_size_manual(name="Variant",values=size_key,drop=TRUE)

variant_linetype_scale_fig2 <- scale_linetype_manual(name="Variant",drop=TRUE,values=linetype_key[c("Original variant","New variant, same kinetics","New variant, different kinetics")])

comparison_color_scale <- scale_color_manual(name="Comparison",drop=TRUE,
                                             values=c("V2-V1,\nsame kinetics"="#009E73",
                                                      "V2-V1,\ndifferent kinetics"="#D55E00",
                                                      "V2-V2,\ndifferent kinetics"="black"))
comparison_fill_scale <- scale_fill_manual(name="Comparison",drop=TRUE,
                                             values=c("V2-V1,\nsame kinetics"="#009E73",
                                                      "V2-V1,\ndifferent kinetics"="#D55E00",
                                                      "V2-V2,\ndifferent kinetics"="black"))
comparison_linetype_scale <- scale_linetype_manual(name="Comparison",drop=TRUE,
                                             values=c("V2-V1,\nsame kinetics"="solid",
                                                      "V2-V1,\ndifferent kinetics"="dashed",
                                                      "V2-V2,\ndifferent kinetics"="solid"))

plot_medians_and_skew <- function(combined_summaries){
  p_medians <- ggplot(combined_summaries) + 
    geom_line(aes(x=t,y=median,col=virus,linetype=virus),size=0.75) +
    variant_color_scale + variant_fill_scale + variant_linetype_scale + 
    scale_y_continuous(trans="reverse") +
    scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
    ylab("Median Ct") +
    theme_overall +
    xlab("Time") +
    theme(legend.position="none") +theme_nice_axes# + theme_no_x_axis
  
  p_skews <- ggplot(combined_summaries) + 
    geom_line(aes(x=t,y=skewness,col=virus,linetype=virus),size=0.75)+
    variant_color_scale + variant_fill_scale + variant_linetype_scale + 
    scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
    ylab("Skewness of Cts") +
    theme_overall +
    xlab("Time") +
    theme(legend.position="none") + theme_nice_axes
  
  return(list(p_medians,p_skews))
}


plot_smooth_mean_cts_symp <- function(ct_values){
  p <- ggplot(ct_values %>% filter(ct < 40)) + 
    geom_smooth(aes(x=sampled_time,y=ct,col=virus,linetype=virus,size=virus,fill=virus)) +
    variant_color_scale + variant_fill_scale + variant_linetype_scale + variant_size_scale +
    scale_y_continuous(trans="reverse") +
    scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
    ylab("Smoothed mean Ct") +
    theme_overall +
    xlab("Time") +
    theme(legend.position="none") +theme_nice_axes
    theme(legend.position="none") + theme_nice_axes
  
  return(p)
}

plot_growth_rate_lineups <- function(ct_combined_summaries){
  p_medians <- ct_combined_summaries %>% pivot_longer(-c(t,virus,gr,inc)) %>% filter(name %in% c("median")) %>% filter(virus != "Overall") %>%
    ggplot() +
    geom_line(aes(x=gr,y=value,col=virus),size=0.75) +
    ylab("Median Ct") +
    xlab("Growth rate") +
    scale_x_continuous(trans="reverse") +
    scale_y_continuous(trans="reverse") +
    variant_color_scale_fig2 + #variant_linetype_scale + 
    theme_overall + theme_nice_axes + theme(legend.position="none")
  
  p_skews <- ct_combined_summaries %>% pivot_longer(-c(t,virus,gr,inc)) %>% filter(name %in% c("skewness")) %>% filter(virus != "Overall") %>%
    ggplot() +
    geom_line(aes(x=gr,y=value,col=virus,size=virus)) +
    scale_x_continuous(trans="reverse") +
    ylab("Skewness of Cts") +
    xlab("Growth rate") +
    variant_color_scale_fig2 +
    theme_overall + #variant_linetype_scale +
    theme_nice_axes + theme(legend.position="bottom")
  return(list(p_medians, p_skews))
}

plot_simulated_ct_curve <- function(vl_pars,ages=1:35,N=1000){
  ## Modal Ct values
  modal_ct <- viral_load_func(vl_pars, ages)
  ## Simulate Ct values
  sim_cts <- virosolver::simulate_viral_loads_example(ages,vl_pars,N)
  p <- ggplot(sim_cts %>% filter(ct < 40)) + 
    geom_rect(ymin=42,ymax=vl_pars["intercept"],xmin=-1,xmax=max(ages)+1,fill="grey70",alpha=0.25) +
    geom_jitter(aes(x=age,y=ct),alpha=0.25,height=0,width=0.25,size=0.2) + 
    geom_line(data=tibble(x=ages,y=modal_ct),aes(x=x,y=y),col=color_key["Original variant"]) +
    geom_line(data=tibble(x=ages,y=modal_ct),aes(x=x,y=y),col=color_key["New variant"],linetype="dotted") +
    coord_flip()+
    coord_cartesian(ylim=c(40,5),xlim=c(0,max(ages))) +
    scale_x_continuous(breaks=seq(0, max(ages),by=5)) +
    xlab("Time since infection") + ylab("Ct value") +
    theme_overall + theme_nice_axes
  return(p)
}

plot_simulated_ct_curve_symptomatic <- function(vl_pars,age_max,N=1000,xmax=20){
  ## Modal Ct values
  ages <- 1:age_max
  ## Simulate Ct values
  sim_cts <- virosolver::simulate_viral_loads_example_symptoms(1:100,vl_pars,N=N)
  p <- ggplot(sim_cts %>% filter(ct_obs < 40)) + 
    geom_rect(ymin=42,ymax=vl_pars["intercept"],xmin=-1,xmax=age_max+1,fill="grey70",alpha=0.25) +
    geom_jitter(aes(x=sampling_delay,y=ct_obs),alpha=0.25,height=0,width=0.25,size=0.2) + 
    geom_smooth(aes(x=sampling_delay,y=ct_obs),
                alpha=0.25,size=0.5) +
    coord_flip()+
    coord_cartesian(ylim=c(40,5),xlim=c(0,xmax)) +
    scale_x_continuous(breaks=seq(0, xmax,by=5)) +
    xlab("Days since symptom onset") + ylab("Ct value") +
    theme_overall + theme_nice_axes
  return(p)
}


plot_indiv_simulated_ct_curve <- function(ct_values,ages=1:35,N=1000){
  ct_values_mean <- ct_values %>% 
    filter(days_since_infection %in% ages) %>%
    filter(ct < 40) %>%
    group_by(days_since_infection,virus) %>% 
    summarize(median_ct=median(ct),mean_ct=mean(ct),mode_ct=calc_mode(ct))
  
  ct_values_tmp <- ct_values %>% 
    filter(days_since_infection %in% ages) %>%
    group_by(virus) %>%
    filter(i %in% sample(unique(i),N))
  
  p <- ggplot(ct_values_tmp) + 
    geom_rect(ymin=42,ymax=vl_pars["intercept"],xmin=-1,xmax=max(ages)+1,fill="grey70",alpha=0.25) +
    geom_jitter(aes(x=days_since_infection,y=ct,col=virus),alpha=0.25,height=0,width=0.25,size=0.2) + 
    geom_line(data=ct_values_mean,aes(x=days_since_infection,y=mode_ct,col=virus,linetype=virus,size=virus)) +
    coord_flip()+
    coord_cartesian(ylim=c(40,5),xlim=c(0,max(ages))) +
    scale_x_continuous(breaks=seq(0, max(ages),by=5)) +
    variant_color_scale_fig2 +
    scale_linetype_manual(name="Variant",values=c("Original variant"="solid","New variant, same kinetics"="dashed", "New variant, different kinetics"="solid")) +
    scale_size_manual(name="Variant",values=c("Original variant"=0.75,"New variant, same kinetics"=1, "New variant, different kinetics"=0.75)) +
    xlab("Time since infection") + ylab("Ct value") +
    theme_overall + theme_nice_axes + theme(legend.position=c(0.7,0.8))
  
  return(p)
}

plot_simulated_ct_curve_symptomatic_OLD <- function(ct_values,age_max=20,N=1000){
  p <- ct_values %>% 
    filter(ct < 40) %>%
    filter(days_since_onset<=age_max) %>% 
    group_by(virus) %>% 
    sample_n(N) %>% 
    ggplot() + 
    geom_rect(ymin=42,ymax=vl_pars["intercept"],xmin=-1,xmax=max(ages)+1,fill="grey70",alpha=0.25) +
    geom_jitter(aes(x=days_since_onset,y=ct,col=virus),height=0,width=0.25,size=0.2,alpha=0.5) +
    geom_smooth(data=ct_dist_symptomatic_dat %>% 
                  filter(ct < 40), aes(x=days_since_onset,y=ct,col=virus,fill=virus),
                alpha=0.25,size=0.5) +
    variant_color_scale + variant_fill_scale +
    coord_flip()+
    coord_cartesian(ylim=c(40,5),xlim=c(0,age_max)) +
    xlab("Days since symptom onset") +
    ylab("Ct value") +
    theme_overall + theme_nice_axes +
    scale_linetype_manual(name="Variant",values=c("Original variant"="solid","New variant, same kinetics"="dashed", 
                                                  "New variant, different kinetics"="solid")) +
    theme(legend.position=c(0.7,0.8))
  return(p)
}




plot_simulated_ct_curve_2variants <- function(vl_pars1, vl_pars2,ages=1:35,N=1000){
  ## Modal Ct values
  modal_ct1 <- tibble(ct=viral_load_func(vl_pars1, ages),age=ages,virus="Original variant")
  modal_ct2 <- tibble(ct=viral_load_func(vl_pars1, ages),age=ages,virus="New variant, same kinetics")
  modal_ct2_alt <- tibble(ct=viral_load_func(vl_pars2, ages),age=ages,virus="New variant, different kinetics")
  modal_ct_all <- bind_rows(modal_ct1,modal_ct2,modal_ct2_alt) %>% mutate(virus=factor(virus,levels=variant_levels))
  
  ## Simulate Ct values
  sim_cts1 <- virosolver::simulate_viral_loads_example(ages,vl_pars1,N) %>% mutate(virus="Original variant")
  sim_cts2 <- virosolver::simulate_viral_loads_example(ages,vl_pars1,N)%>% mutate(virus="New variant, same kinetics")
  sim_cts2_alt <- virosolver::simulate_viral_loads_example(ages,vl_pars2,N)%>% mutate(virus="New variant, different kinetics")
  sim_cts_all <- bind_rows(sim_cts1, sim_cts2, sim_cts2_alt) %>% mutate(virus=factor(virus,levels=variant_levels))
  
  p <- ggplot(sim_cts_all %>% filter(ct < 40)) + 
    geom_rect(ymin=42,ymax=vl_pars["intercept"],xmin=-1,xmax=max(ages)+1,fill="grey70",alpha=0.25) +
    geom_jitter(aes(x=age,y=ct,col=virus),alpha=0.25,height=0,width=0.25,size=0.2) + 
    geom_line(data=modal_ct_all,aes(x=age,y=ct,col=virus,linetype=virus)) +
    coord_flip()+
    coord_cartesian(ylim=c(40,5),xlim=c(0,max(ages))) +
    scale_x_continuous(breaks=seq(0, max(ages),by=5)) +
    variant_color_scale_fig2 +
    scale_linetype_manual(name="Variant",values=c("Original variant"="solid","New variant, same kinetics"="dashed", "New variant, different kinetics"="solid")) +
    scale_size_manual(name="Variant",values=c("Original variant"=0.75,"New variant, same kinetics"=1, "New variant, different kinetics"=0.75)) +
    xlab("Time since infection") + ylab("Ct value") +
    theme_overall + theme_nice_axes + theme(legend.position=c(0.7,0.8))
  return(p)
}


plot_simulated_ct_curve_2variants_symptomatic <- function(vl_pars1, vl_pars2,N=1000,xmax=20){
  ## Modal Ct values
  ages <- 1:100
  
  ## Modal Ct values
  modal_ct1 <- tibble(ct=viral_load_func(vl_pars1, ages),age=ages,virus="Original variant")
  modal_ct2 <- tibble(ct=viral_load_func(vl_pars1, ages),age=ages,virus="New variant, same kinetics")
  modal_ct2_alt <- tibble(ct=viral_load_func(vl_pars2, ages),age=ages,virus="New variant, different kinetics")
  modal_ct_all <- bind_rows(modal_ct1,modal_ct2,modal_ct2_alt) %>% mutate(virus=factor(virus,levels=variant_levels))
  
  ## Simulate Ct values
  sim_cts1 <- virosolver::simulate_viral_loads_example_symptoms(ages,vl_pars1,N=N) %>% mutate(virus="Original variant")
  sim_cts2 <- virosolver::simulate_viral_loads_example_symptoms(ages,vl_pars1,N=N)%>% mutate(virus="New variant, same kinetics")
  sim_cts2_alt <- virosolver::simulate_viral_loads_example_symptoms(ages,vl_pars2,N=N)%>% mutate(virus="New variant, different kinetics")
  sim_cts_all <- bind_rows(sim_cts1, sim_cts2, sim_cts2_alt) %>% mutate(virus=factor(virus,levels=variant_levels))
  
  p <- ggplot(sim_cts_all %>% filter(ct_obs < 40)) + 
    geom_rect(ymin=42,ymax=vl_pars["intercept"],xmin=-1,xmax=max(ages)+1,fill="grey70",alpha=0.25) +
    geom_jitter(aes(x=sampling_delay,y=ct_obs,col=virus),alpha=0.25,height=0,width=0.25,size=0.2) + 
    geom_smooth(aes(x=sampling_delay,y=ct_obs,col=virus,linetype=virus,fill=virus),alpha=0.1,size=0.5) +
    coord_flip()+
    coord_cartesian(ylim=c(40,5),xlim=c(0,xmax)) +
    scale_x_continuous(breaks=seq(0, xmax,by=5)) +
    variant_color_scale_fig2 + variant_fill_scale_fig2 +
    scale_linetype_manual(name="Variant",values=c("Original variant"="solid","New variant, same kinetics"="dashed", "New variant, different kinetics"="solid")) +
    scale_size_manual(name="Variant",values=c("Original variant"=0.75,"New variant, same kinetics"=1, "New variant, different kinetics"=0.75)) +
    xlab("Days since symptom onset") + ylab("Ct value") +
    theme_overall + theme_nice_axes + theme(legend.position=c(0.7,0.8))
  return(p)
}



p_sim_ct_compare_naive <- function(vl_pars1,vl_pars2,virus1_inc,virus2_inc, ages,samp_time,N=100,dotsize=1,symptom_surveillance=FALSE) {
  cts_1 <- tibble(ct=simulate_cross_section(vl_pars1, ages, virus1_inc,obs_time=samp_time,N=N,use_pos=TRUE,symptom_surveillance=symptom_surveillance),virus="Original variant")
  cts_2 <- tibble(ct=simulate_cross_section(vl_pars1, ages, virus2_inc,obs_time=samp_time,N=N,use_pos=TRUE,symptom_surveillance=symptom_surveillance),virus="New variant, same kinetics")
  cts_2_alt <- tibble(ct=simulate_cross_section(vl_pars2, ages, virus2_inc,obs_time=samp_time,N=N,use_pos=TRUE,symptom_surveillance=symptom_surveillance),virus="New variant, different kinetics")
  cts_sim_comb <- bind_rows(cts_1,cts_2,cts_2_alt) %>% mutate(virus=factor(virus,levels=variant_levels))
  
  pval1 <- as.numeric(wilcox.test(cts_1$ct,cts_2$ct,alternative="two.sided")["p.value"])
  pval3 <- as.numeric(wilcox.test(cts_1$ct,cts_2_alt$ct,alternative="two.sided")["p.value"])
  pval2 <- as.numeric(wilcox.test(cts_2$ct,cts_2_alt$ct,alternative="two.sided")["p.value"])
  
  pval_labels <- tibble(x=c(1.5,2.5,2),y=c(12,10,8)-1,pval=signif(c(pval1,pval2,pval3),3),
                        signf=ifelse(pval>0.05,"",ifelse(pval > 0.01,"*",ifelse(pval>0.001,"**","***"))),
                        label=paste0(signf,"p=",pval))
  pval_lines <- tibble(x=c(1,2,1),xend=c(2,3,3),y=c(12,10,8),yend=c(12,10,8))
  p_ct_samp <- ggplot(cts_sim_comb) + 
    geom_violin(aes(x=virus,y=ct,fill=virus),alpha=0.1,trim=FALSE,col=NA,width=0.75) +
    geom_dotplot(aes(x=virus,y=ct,fill=virus),binaxis="y",stackdir="center",binwidth=1,dotsize=dotsize) + 
    geom_errorbar(data=cts_sim_comb%>%group_by(virus)%>%summarize(median_ct=median(ct)),
                  aes(y=median_ct,ymin=median_ct,ymax=median_ct,x=virus),size=0.5,col="grey10") +
    geom_segment(data=pval_lines,aes(x=x,xend=xend,y=y,yend=yend),size=0.5) +
    geom_text(data=pval_labels,aes(x=x,y=y,label=label),size=2) +
    variant_fill_scale + variant_color_scale +
    scale_y_continuous(trans="reverse",limits=c(40,7)) +
    ylab("Ct value") +
    theme_overall + theme_nice_axes + theme(legend.position="none",axis.title.x=element_blank(),plot.title=element_text(size=7)) +
    scale_x_discrete(labels = c("Original variant" = "Original variant",
                                "New variant, same kinetics"="New variant,\nsame kinetics",
                                "New variant, different kinetics"="New variant,\ndifferent kinetics")) +
    ggtitle(paste0("Wilcoxon Rank Sum test on samples from day ", samp_time))
  p_ct_samp
}


p_sim_ct_compare_growth <- function(vl_pars1,vl_pars2,virus1_inc,virus2_inc, ages,combined_summaries,growth_rate_samp,N=100,dotsize=1,symptom_surveillance=FALSE) {
  samp_times <- ct_combined_summaries %>% 
    group_by(virus) %>% 
    mutate(gr_diff=abs(gr - growth_rate_samp)) %>% 
    filter(gr_diff==min(gr_diff))
  
  samp_time1 <- samp_times %>% filter(virus == "Original variant") %>% pull(t)
  samp_time2 <- samp_times %>% filter(virus == "New variant, same kinetics") %>% pull(t)
  samp_time2_alt <- samp_times %>% filter(virus == "New variant, different kinetics") %>% pull(t)
  
  cts_1 <- tibble(ct=simulate_cross_section(vl_pars1, ages, virus1_inc,obs_time=samp_time1,N=N,use_pos=TRUE,symptom_surveillance=symptom_surveillance),virus="Original variant")
  cts_2 <- tibble(ct=simulate_cross_section(vl_pars1, ages, virus2_inc,obs_time=samp_time2,N=N,use_pos=TRUE,symptom_surveillance=symptom_surveillance),virus="New variant, same kinetics")
  cts_2_alt <- tibble(ct=simulate_cross_section(vl_pars2, ages, virus2_inc,obs_time=samp_time2_alt,N=N,use_pos=TRUE,symptom_surveillance=symptom_surveillance),virus="New variant, different kinetics")
  cts_sim_comb <- bind_rows(cts_1,cts_2,cts_2_alt) %>% mutate(virus=factor(virus,levels=variant_levels))
  
  pval1 <- as.numeric(wilcox.test(cts_1$ct,cts_2$ct,alternative="two.sided")["p.value"])
  pval3 <- as.numeric(wilcox.test(cts_1$ct,cts_2_alt$ct,alternative="two.sided")["p.value"])
  pval2 <- as.numeric(wilcox.test(cts_2$ct,cts_2_alt$ct,alternative="two.sided")["p.value"])
  
  pval_labels <- tibble(x=c(1.5,2.5,2),y=c(12,10,8)-1,pval=signif(c(pval1,pval2,pval3),3),
                        signf=ifelse(pval>0.05,"",ifelse(pval > 0.01,"*",ifelse(pval>0.001,"**","***"))),
                        label=paste0(signf,"p=",pval))
  pval_lines <- tibble(x=c(1,2,1),xend=c(2,3,3),y=c(12,10,8),yend=c(12,10,8))
  p_ct_samp <- ggplot(cts_sim_comb) + 
    geom_violin(aes(x=virus,y=ct,fill=virus),alpha=0.1,trim=FALSE,col=NA,width=0.5) +
    geom_dotplot(aes(x=virus,y=ct,fill=virus),binaxis="y",stackdir="center",binwidth=1,dotsize=dotsize) + 
    geom_errorbar(data=cts_sim_comb%>%group_by(virus)%>%summarize(median_ct=median(ct)),
                  aes(y=median_ct,ymin=median_ct,ymax=median_ct,x=virus),size=0.5,col="grey10") +
    geom_segment(data=pval_lines,aes(x=x,xend=xend,y=y,yend=yend),size=0.5) +
    geom_text(data=pval_labels,aes(x=x,y=y,label=label),size=2) +
    variant_fill_scale + variant_color_scale +
    scale_y_continuous(trans="reverse",limits=c(40,7)) +
    ylab("Ct value") +
    theme_overall + theme_nice_axes + theme(legend.position="none",axis.title.x=element_blank(),plot.title=element_text(size=7)) +
    scale_x_discrete(labels = c("Original variant" = "Original variant",
                                "New variant, same kinetics"="New variant,\nsame kinetics",
                                "New variant, different kinetics"="New variant,\ndifferent kinetics")) +
    ggtitle(paste0("Wilcoxon Rank Sum test on samples from growth rate=", growth_rate_samp))
  p_ct_samp
}

internal_plot_power_compare_power_symp <- function(all_results, all_results_summary, true_peak_diff, samp_sizes,ver="wilcoxon"){
  p_differences <- ggplot(all_results %>% mutate(samp_size_label=factor(samp_size_label,levels=paste0("N=",samp_sizes)))) + 
    #geom_jitter(aes(x=comparison,y=difference,col=comparison),height=0,width=0.25,size=0.25,alpha=0.5) + 
    geom_boxplot(aes(x=comparison,y=difference,col=comparison,fill=comparison),size=0.25,alpha=0.25) +
    geom_hline(yintercept=0,linetype="dashed",col="#009E73") + 
    geom_hline(yintercept=true_peak_diff,linetype="dashed",col="#D55E00") + 
    comparison_color_scale + comparison_fill_scale + 
    facet_wrap(~samp_size_label,nrow=1)+
    ylab("Median Ct difference") +
    xlab("") +
    theme_overall + 
    theme_nice_axes +
    theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="none") +
    labs(tag="A")
  
  if(ver == "regression"){
    p_differences <- p_differences + ylab("Coefficient for virus variable")
  }
  
  p_t1_error <- ggplot(all_results_summary %>% 
                         filter(comparison %in% c("V2-V1,\ndifferent kinetics"))%>%
                         mutate(samp_size=factor(samp_size,levels=samp_sizes))) + 
    geom_hline(yintercept=0.95,linetype="dashed",alpha=0.5) +
    geom_line(aes(x=as.numeric(samp_size),y=prop_different,col=comparison)) + 
    #comparison_color_scale + comparison_linetype_scale +
    scale_color_manual(name="Comparison",values=c("V2-V1,\ndifferent kinetics"="#D55E00")) +
    ylab("Power") +
    xlab("") +
    theme_overall + 
    theme_nice_axes +
    scale_y_continuous(limits=c(0,1)) +
    theme(legend.position=c(0.8,0.4)) +
    scale_x_continuous(labels=samp_sizes,breaks=seq_along(samp_sizes)) +
    labs(tag="B")
  
  p_t2_error <- ggplot(all_results_summary %>% filter(comparison %in% c("V2-V1,\nsame kinetics"))%>%
                         mutate(samp_size=factor(samp_size,levels=samp_sizes))) + 
    geom_hline(yintercept=0.05,linetype="dashed",alpha=0.5) +
    geom_line(aes(x=as.numeric(samp_size),y=prop_different,col=comparison)) + 
    comparison_color_scale +
    scale_color_manual(name="Comparison",values=c("V2-V1,\nsame kinetics"="#009E73")) +
    scale_y_continuous(limits=c(0,1)) + 
    ylab("Type 1 error") +
    xlab("Sample size") +
    theme_overall + 
    theme_nice_axes +
    theme(legend.position=c(0.4,0.8)) +
    scale_x_continuous(labels=samp_sizes,breaks=seq_along(samp_sizes)) +
    labs(tag="C")
  
  p_main <- p_differences/p_t1_error/p_t2_error 
  p_main
}


p_sim_ct_compare_power <- function(vl_pars1,vl_pars2,virus1_inc,virus2_inc, 
                                   ages,samp_time,trials=100,samp_sizes=100,
                                   alpha=0.05,
                                   align_gr=FALSE,grs=NULL,gr_test=NULL) {
  if(align_gr){
    samp_times <- grs %>% 
      group_by(virus) %>% 
      mutate(gr_diff=abs(gr - gr_test)) %>% 
      filter(gr_diff==min(gr_diff,na.rm=TRUE))
    
    samp_time1 <- samp_times %>% filter(virus == "Original variant") %>% pull(t)
    samp_time2 <- samp_times %>% filter(virus == "New variant, same kinetics") %>% pull(t)
    samp_time2_alt <- samp_times %>% filter(virus == "New variant, different kinetics") %>% pull(t)
  } else {
    samp_time1 <- samp_time2 <- samp_time2_alt <- samp_time
  }
  
  all_results <- list()
  for(i in seq_along(samp_sizes)){
    N <- samp_sizes[i]
    print(paste0("Sample size: ", N))
    cts_1 <- simulate_m_cross_sections(vl_pars1, ages, virus1_inc,obs_time=samp_time1,N=N,m=trials) %>% mutate(virus="Original variant")
    cts_2 <- simulate_m_cross_sections(vl_pars1, ages, virus2_inc,obs_time=samp_time2,N=N,m=trials) %>% mutate(virus="New variant, same kinetics")
    cts_2_alt <- simulate_m_cross_sections(vl_pars2, ages, virus2_inc,obs_time=samp_time2_alt,N=N,m=trials) %>% mutate(virus="New variant, different kinetics")
    cts_sim_comb <- bind_rows(cts_1,cts_2,cts_2_alt) %>% mutate(virus=factor(virus,levels=variant_levels))
    
    true_peak_diff <- vl_pars2["viral_peak"] - vl_pars1["viral_peak"]
    true_overall_diff_1v2 <- median(cts_2$ct) - median(cts_1$ct)
    true_overall_diff_1v2b <- median(cts_2_alt$ct) - median(cts_1$ct)
    true_overall_diff_2v2b <- median(cts_2_alt$ct) - median(cts_2$ct)
    
    results <- list()
    for(trial in unique(cts_sim_comb$sampno)){
      v1_cts <- cts_sim_comb %>% filter(sampno == trial,virus=="Original variant") %>% pull(ct)
      v2_cts <- cts_sim_comb %>% filter(sampno == trial,virus=="New variant, same kinetics") %>% pull(ct)
      v2_alt_cts <- cts_sim_comb %>% filter(sampno == trial,virus=="New variant, different kinetics") %>% pull(ct)
      
      diff1v2 <- median(v2_cts) - median(v1_cts)
      diff1v2alt <- median(v2_alt_cts) - median(v1_cts)
      diff2v2alt <- median(v2_alt_cts) - median(v2_cts)
      
      test1v2 <- as.numeric(wilcox.test(v1_cts,v2_cts,alternative="two.sided")["p.value"])
      test1v2alt <- as.numeric(wilcox.test(v1_cts,v2_alt_cts,alternative="two.sided")["p.value"])
      test2v2alt <- as.numeric(wilcox.test(v2_alt_cts,v2_cts,alternative="two.sided")["p.value"])
      
      results[[trial]] <- tibble(comparison=c("V2-V1,\nsame kinetics","V2-V1,\ndifferent kinetics","V2-V2,\ndifferent kinetics"),
                                 difference=c(diff1v2,diff1v2alt,diff2v2alt),
                                 pvalue=c(test1v2,test1v2alt,test2v2alt),trial=trial) %>% 
        mutate(signif=as.numeric(pvalue<alpha),samp_size=N,samp_size_label=paste0("N=",N))
      
    }
    results <- do.call("bind_rows",results)
    all_results[[N]] <- results
  }
  all_results <- do.call("bind_rows",all_results)
  all_results_summary <- all_results %>% group_by(comparison,samp_size) %>% summarize(prop_different=sum(signif)/n())
  
  all_results$comparison <- factor(all_results$comparison, levels=c("V2-V1,\nsame kinetics","V2-V1,\ndifferent kinetics","V2-V2,\ndifferent kinetics"))
  all_results_summary$comparison <- factor(all_results_summary$comparison, levels=c("V2-V1,\nsame kinetics","V2-V1,\ndifferent kinetics","V2-V2,\ndifferent kinetics"))
  
  p_main <- internal_plot_power_compare_power_symp(all_results, all_results_summary, true_peak_diff, samp_sizes)
  
  list(all_results,all_results_summary,p_main)
}


p_sim_ct_indiv_compare_power <- function(cts_indiv_comb,samp_time,trials=100,samp_sizes=100,alpha=0.05,
                                   align_gr=FALSE,grs=NULL,gr_test=NULL,true_peak_diff=0) {
 
  
  if(align_gr){
    samp_times <- grs %>% 
      group_by(virus) %>% 
      mutate(gr_diff=abs(gr - gr_test)) %>% 
      filter(gr_diff==min(gr_diff,na.rm=TRUE))
    
    samp_time1 <- samp_times %>% filter(virus == "Original variant") %>% pull(t)
    samp_time2 <- samp_times %>% filter(virus == "New variant, same kinetics") %>% pull(t)
    samp_time2_alt <- samp_times %>% filter(virus == "New variant, different kinetics") %>% pull(t)
  } else {
    samp_time1 <- samp_time2 <- samp_time2_alt <- samp_time
  }
  all_results <- list()
  for(i in seq_along(samp_sizes)){
    N <- samp_sizes[i]
    print(paste0("Sample size: ", N))
    ## Generate fake Ct distributions for the desired time and virus
    ## Original variant
    cts_1_tmp <- cts_indiv_comb %>% filter(t == samp_time1, virus=="Original variant") %>% pull(ct)
    cts_1_all <- resample_ct_dist(cts_1_tmp,samp_size=N,N=trials,with_replacement=TRUE,bootstrap_cts=FALSE,cts=seq(0,40,by=0.1))
    ## New variant, same kinetics
    cts_2_tmp <- cts_indiv_comb %>% filter(t == samp_time2, virus=="New variant, same kinetics") %>% pull(ct)
    cts_2_all <- resample_ct_dist(cts_2_tmp,samp_size=N,N=trials,with_replacement=TRUE,bootstrap_cts=FALSE,cts=seq(0,40,by=0.1))
    ## Original variant
    cts_2_alt_tmp <- cts_indiv_comb %>% filter(t == samp_time2_alt, virus=="New variant, different kinetics") %>% pull(ct)
    cts_2_alt_all <- resample_ct_dist(cts_2_alt_tmp,samp_size=N,N=trials,with_replacement=TRUE,bootstrap_cts=FALSE,cts=seq(0,40,by=0.1))
    
    results <- list()
    for(j in 1:trials){
      v1_cts <- cts_1_all %>% filter(trial == j) %>% pull(ct)
      v2_cts <- cts_2_all %>% filter(trial == j) %>% pull(ct)
      v2_alt_cts <- cts_2_alt_all %>% filter(trial == j) %>% pull(ct)
      
      diff1v2 <- median(v2_cts) - median(v1_cts)
      diff1v2alt <- median(v2_alt_cts) - median(v1_cts)
      diff2v2alt <- median(v2_alt_cts) - median(v2_cts)
      
      test1v2 <- as.numeric(wilcox.test(v1_cts,v2_cts,alternative="two.sided")["p.value"])
      test1v2alt <- as.numeric(wilcox.test(v1_cts,v2_alt_cts,alternative="two.sided")["p.value"])
      test2v2alt <- as.numeric(wilcox.test(v2_alt_cts,v2_cts,alternative="two.sided")["p.value"])
      
      results[[j]] <- tibble(comparison=c("V2-V1,\nsame kinetics","V2-V1,\ndifferent kinetics","V2-V2,\ndifferent kinetics"),
                                 difference=c(diff1v2,diff1v2alt,diff2v2alt),
                                 pvalue=c(test1v2,test1v2alt,test2v2alt),trial=j) %>% 
        mutate(signif=as.numeric(pvalue<alpha),samp_size=N,samp_size_label=paste0("N=",N))
    }
    results <- do.call("bind_rows",results)
    all_results[[N]] <- results
  }
  all_results <- do.call("bind_rows",all_results)
  all_results_summary <- all_results %>% group_by(comparison,samp_size) %>% summarize(prop_different=sum(signif)/n())
  
  all_results$comparison <- factor(all_results$comparison, levels=c("V2-V1,\nsame kinetics","V2-V1,\ndifferent kinetics","V2-V2,\ndifferent kinetics"))
  all_results_summary$comparison <- factor(all_results_summary$comparison, levels=c("V2-V1,\nsame kinetics","V2-V1,\ndifferent kinetics","V2-V2,\ndifferent kinetics"))
  
  p_main <- internal_plot_power_compare_power_symp(all_results, all_results_summary, true_peak_diff, samp_sizes)
  
  list(all_results,all_results_summary,p_main)
}


p_sim_ct_compare_power_symp <- function(ct_dist_symptomatic_dat, samp_time,trials=100,
                                        samp_sizes=100,alpha=0.05,true_peak_diff=0,samp_window=7) {
  samp_time1 <- samp_time2 <- samp_time2_alt <- samp_time
  
  all_results <- list()
  for(i in seq_along(samp_sizes)){
    N <- samp_sizes[i]
    print(paste0("Sample size: ", N))
    ## Generate fake Ct distributions for the desired time and virus
    ## Original variant
    cts_1_tmp <- ct_dist_symptomatic_dat %>% filter(sampled_time >= samp_time1-(samp_window/2),sampled_time <= samp_time1+(samp_window/2), 
                                                    virus=="Original variant") %>% pull(ct)
    cts_1_all <- resample_ct_dist(cts_1_tmp,samp_size=N,N=trials,with_replacement=TRUE,bootstrap_cts=FALSE,cts=seq(0,40,by=0.1))
    ## New variant, same kinetics
    cts_2_tmp <- ct_dist_symptomatic_dat %>% filter(sampled_time >= samp_time2-(samp_window/2),sampled_time <= samp_time2+(samp_window/2), 
                                                    virus=="New variant, same kinetics") %>% pull(ct)
    cts_2_all <- resample_ct_dist(cts_2_tmp,samp_size=N,N=trials,with_replacement=TRUE,bootstrap_cts=FALSE,cts=seq(0,40,by=0.1))
    ## Original variant
    cts_2_alt_tmp <- ct_dist_symptomatic_dat %>% filter(sampled_time >= samp_time2_alt-(samp_window/2),sampled_time <= samp_time2_alt+(samp_window/2), 
                                                        virus=="New variant, different kinetics") %>% pull(ct)
    cts_2_alt_all <- resample_ct_dist(cts_2_alt_tmp,samp_size=N,N=trials,with_replacement=TRUE,bootstrap_cts=FALSE,cts=seq(0,40,by=0.1))
    
    results <- list()
    for(j in 1:trials){
      v1_cts <- cts_1_all %>% filter(trial == j) %>% pull(ct)
      v2_cts <- cts_2_all %>% filter(trial == j) %>% pull(ct)
      v2_alt_cts <- cts_2_alt_all %>% filter(trial == j) %>% pull(ct)
      
      diff1v2 <- median(v2_cts) - median(v1_cts)
      diff1v2alt <- median(v2_alt_cts) - median(v1_cts)
      diff2v2alt <- median(v2_alt_cts) - median(v2_cts)
      
      test1v2 <- as.numeric(wilcox.test(v1_cts,v2_cts,alternative="two.sided")["p.value"])
      test1v2alt <- as.numeric(wilcox.test(v1_cts,v2_alt_cts,alternative="two.sided")["p.value"])
      test2v2alt <- as.numeric(wilcox.test(v2_alt_cts,v2_cts,alternative="two.sided")["p.value"])
      
      results[[j]] <- tibble(comparison=c("V2-V1,\nsame kinetics","V2-V1,\ndifferent kinetics","V2-V2,\ndifferent kinetics"),
                             difference=c(diff1v2,diff1v2alt,diff2v2alt),
                             pvalue=c(test1v2,test1v2alt,test2v2alt),trial=j) %>% 
        mutate(signif=as.numeric(pvalue<alpha),samp_size=N,samp_size_label=paste0("N=",N))
    }
    results <- do.call("bind_rows",results)
    all_results[[N]] <- results
  }
  all_results <- do.call("bind_rows",all_results)
  all_results_summary <- all_results %>% group_by(comparison,samp_size) %>% summarize(prop_different=sum(signif)/n())
  
  all_results$comparison <- factor(all_results$comparison, levels=c("V2-V1,\nsame kinetics","V2-V1,\ndifferent kinetics","V2-V2,\ndifferent kinetics"))
  all_results_summary$comparison <- factor(all_results_summary$comparison, levels=c("V2-V1,\nsame kinetics","V2-V1,\ndifferent kinetics","V2-V2,\ndifferent kinetics"))
  
  p_main <- internal_plot_power_compare_power_symp(all_results, all_results_summary, true_peak_diff, samp_sizes)
  
  list(all_results,all_results_summary,p_main)
}



p_sim_ct_compare_power_symp_regression <- function(ct_dist_symptomatic_dat, samp_time,trials=100,
                                                   samp_sizes=100,alpha=0.05,true_peak_diff=0,samp_window=7) {
  samp_time1 <- samp_time2 <- samp_time2_alt <- samp_time
  
  all_results <- list()
  all_results_wilcox <- list()
  
  for(i in seq_along(samp_sizes)){
    N <- samp_sizes[i]
    print(paste0("Sample size: ", N))
    
    results <- list()
    results_wilcox <- list()
    
    for(j in 1:trials){
      ## Resample Ct values from the time range for each variant
      v1_cts <- ct_dist_symptomatic_dat %>% filter(sampled_time >= samp_time1-(samp_window/2),sampled_time <= samp_time1+(samp_window/2),
                                                   virus=="Original variant") %>% sample_n(N,replace=TRUE)
      v2_cts <- ct_dist_symptomatic_dat %>% filter(sampled_time >= samp_time2-(samp_window/2),sampled_time <= samp_time2+(samp_window/2),
                                                   virus=="New variant, same kinetics") %>% sample_n(N,replace=TRUE)
      v2_alt_cts <- ct_dist_symptomatic_dat %>% filter(sampled_time >= samp_time2_alt-(samp_window/2),sampled_time <= samp_time2_alt+(samp_window/2),
                                                       virus=="New variant, different kinetics") %>% sample_n(N,replace=TRUE)
      
      
      ## Pull Ct values and find difference in medians, Wilcoxon test
      v1_cts_wilcox <- v1_cts %>% pull(ct)
      v2_cts_wilcox <- v2_cts %>% pull(ct)
      v2_alt_cts_wilcox <- v2_alt_cts %>% pull(ct)
      
      diff1v2 <- median(v2_cts_wilcox) - median(v1_cts_wilcox)
      diff1v2alt <- median(v2_alt_cts_wilcox) - median(v1_cts_wilcox)
      diff2v2alt <- median(v2_alt_cts_wilcox) - median(v2_cts_wilcox)
      
      test1v2 <- as.numeric(wilcox.test(v1_cts_wilcox,v2_cts_wilcox,alternative="two.sided")["p.value"])
      test1v2alt <- as.numeric(wilcox.test(v1_cts_wilcox,v2_alt_cts_wilcox,alternative="two.sided")["p.value"])
      test2v2alt <- as.numeric(wilcox.test(v2_alt_cts_wilcox,v2_cts_wilcox,alternative="two.sided")["p.value"])
      
      results_wilcox[[j]] <- tibble(comparison=c("V2-V1,\nsame kinetics","V2-V1,\ndifferent kinetics","V2-V2,\ndifferent kinetics"),
                             difference=c(diff1v2,diff1v2alt,diff2v2alt),
                             pvalue=c(test1v2,test1v2alt,test2v2alt),trial=j) %>% 
        mutate(signif=as.numeric(pvalue<alpha),samp_size=N,samp_size_label=paste0("N=",N))
      
      
      ## See how much variation "variant" describes in linear regression model
      v_comp_1 <- bind_rows(v1_cts, v2_cts)
      fit1 <- lm(ct ~ days_since_onset + virus,data=v_comp_1)
      
      v_comp_2 <- bind_rows(v1_cts, v2_alt_cts)
      fit2 <- lm(ct ~ days_since_onset + virus,data=v_comp_2)
      
      v_comp_3 <- bind_rows(v2_alt_cts, v2_cts)
      fit3 <- lm(ct ~ days_since_onset + virus,data=v_comp_3)
      
      ## Difference due to virus
      diff1v2 <- summary(fit1)$coefficients[3,"Estimate"]
      diff1v2alt <- summary(fit2)$coefficients[3,"Estimate"]
      diff2v2alt <- summary(fit3)$coefficients[3,"Estimate"]
      
      test1v2 <- summary(fit1)$coefficients[3,4]
      test1v2alt <- summary(fit2)$coefficients[3,4]
      test2v2alt <- summary(fit3)$coefficients[3,4]
      
      lower_ci1 <- confint(fit1)[3,1]
      lower_ci2 <- confint(fit2)[3,1]
      lower_ci3 <- confint(fit3)[3,1]
      
      upper_ci1 <- confint(fit1)[3,2]
      upper_ci2 <- confint(fit2)[3,2]
      upper_ci3 <- confint(fit3)[3,2]
      
      
      ## Difference due to days since onset
      onset_diff1v2 <- summary(fit1)$coefficients[2,"Estimate"]
      onset_diff1v2alt <- summary(fit2)$coefficients[2,"Estimate"]
      onset_diff2v2alt <- summary(fit3)$coefficients[2,"Estimate"]
      
      onset_test1v2 <- summary(fit1)$coefficients[2,4]
      onset_test1v2alt <- summary(fit2)$coefficients[2,4]
      onset_test2v2alt <- summary(fit3)$coefficients[2,4]
      
      onset_lower_ci1 <- confint(fit1)[2,1]
      onset_lower_ci2 <- confint(fit2)[2,1]
      onset_lower_ci3 <- confint(fit3)[2,1]
      
      onset_upper_ci1 <- confint(fit1)[2,2]
      onset_upper_ci2 <- confint(fit2)[2,2]
      onset_upper_ci3 <- confint(fit3)[2,2]
      
      results[[j]] <- tibble(comparison=c("V2-V1,\nsame kinetics","V2-V1,\ndifferent kinetics","V2-V2,\ndifferent kinetics"),
                             difference=c(diff1v2,diff1v2alt,diff2v2alt),
                             pvalue=c(test1v2,test1v2alt,test2v2alt),
                             lower_confint=c(lower_ci1,lower_ci2,lower_ci3),
                             upper_confint=c(upper_ci1,upper_ci2,upper_ci3),
                             
                             onset_difference=c(onset_diff1v2,onset_diff1v2alt,onset_diff2v2alt),
                             onset_pvalue=c(onset_test1v2,onset_test1v2alt,onset_test2v2alt),
                             onset_lower_confint=c(onset_lower_ci1,onset_lower_ci2,onset_lower_ci3),
                             onset_upper_confint=c(onset_upper_ci1,onset_upper_ci2,onset_upper_ci3),
                             
                             true_diff=c(0, true_peak_diff, true_peak_diff),
                             correct=c(lower_ci1 < 0 & upper_ci1 > 0, upper_ci2 < 0, upper_ci3 < 0),
                             trial=j) %>% 
        mutate(signif=as.numeric(pvalue<alpha),samp_size=N,samp_size_label=paste0("N=",N))
    }
    results <- do.call("bind_rows",results)
    results_wilcox <- do.call("bind_rows",results_wilcox)
    all_results[[N]] <- results
    all_results_wilcox[[N]] <- results_wilcox
  }
  
  ## Linear regression model results
  all_results <- do.call("bind_rows",all_results)
  all_results_summary <- all_results %>% group_by(comparison,samp_size) %>% summarize(prop_correct=sum(correct)/n()) %>%
    mutate(prop_different=ifelse(comparison == "V2-V1,\nsame kinetics",1-prop_correct,prop_correct))
  
  all_results$comparison <- factor(all_results$comparison, 
                                   levels=c("V2-V1,\nsame kinetics","V2-V1,\ndifferent kinetics","V2-V2,\ndifferent kinetics"))
  all_results_summary$comparison <- factor(all_results_summary$comparison, 
                                           levels=c("V2-V1,\nsame kinetics","V2-V1,\ndifferent kinetics","V2-V2,\ndifferent kinetics"))
  
  ## Wilcoxon test results
  all_results_wilcox <- do.call("bind_rows",all_results_wilcox)
  all_results_summary_wilcox <- all_results_wilcox %>% group_by(comparison,samp_size) %>% summarize(prop_different=sum(signif)/n())
  
  all_results_wilcox$comparison <- factor(all_results_wilcox$comparison, levels=c("V2-V1,\nsame kinetics","V2-V1,\ndifferent kinetics","V2-V2,\ndifferent kinetics"))
  all_results_summary_wilcox$comparison <- factor(all_results_summary_wilcox$comparison, levels=c("V2-V1,\nsame kinetics","V2-V1,\ndifferent kinetics","V2-V2,\ndifferent kinetics"))
  
  p_main_lm <- internal_plot_power_compare_power_symp(all_results, all_results_summary,true_peak_diff, samp_sizes,ver="regression")
  p_main_wilcox <- internal_plot_power_compare_power_symp(all_results_wilcox, all_results_summary_wilcox,true_peak_diff, samp_sizes)
  list(all_results,all_results_summary,
       all_results_wilcox, all_results_summary_wilcox,
       p_main_lm, p_main_wilcox)
}

p_compare_estimated_curves <- function(chain,ages=1:35,N=1000, nsamp=100,true_pars1=NULL,true_pars2=NULL){
  ## If passed to function, get the true modal viral kinetics curve
  if(!is.null(true_pars1)){
    vl1 <- viral_load_func(true_pars1, ages)
    vl2 <- viral_load_func(true_pars2, ages)
    
    ## Modal Ct values
    true_modal_ct1 <- tibble(ct=vl1,age=ages,virus="Original variant")
    true_modal_ct2 <- tibble(ct=vl2,age=ages,virus="New variant")
    true_modal_cts <- bind_rows(true_modal_ct1,true_modal_ct2)
  }
  
  
  samps <- sample(chain$sampno,size=nsamp)
  modal_ct_all <- NULL
  sim_cts_all <- NULL
  for(i in seq_along(samps)){
    samp <- samps[i]
    
    pars_tmp <- get_index_pars(chain, samp)
    pars_tmp_variant <- pars_tmp
    
    ## Viral load control points
    pars_tmp_variant["viral_peak"] <- pars_tmp_variant["viral_peak"]*pars_tmp_variant["viral_peak_scale"]
    pars_tmp_variant["level_switch"] <- pars_tmp_variant["level_switch"]*pars_tmp_variant["level_switch_scale"]
    ## Timing pars
    pars_tmp_variant["desired_mode"] <- pars_tmp_variant["desired_mode"]*pars_tmp_variant["desired_mode_scale"]
    pars_tmp_variant["t_switch"] <- pars_tmp_variant["t_switch"]*pars_tmp_variant["t_switch_scale"]
    ## Detect loss
    pars_tmp_variant["prob_detect"] <- pars_tmp_variant["prob_detect"]*pars_tmp_variant["prob_detect_scale"]
    ## Measurement parameters
    pars_tmp_variant["obs_sd"] <- pars_tmp_variant["obs_sd"]*pars_tmp_variant["obs_sd_scale"]
    pars_tmp_variant["sd_mod"] <- pars_tmp_variant["sd_mod"]*pars_tmp_variant["sd_mod_scale"]
    pars_tmp_variant["sd_mod_wane"] <- pars_tmp_variant["sd_mod_wane"]*pars_tmp_variant["sd_mod_wane_scale"]
    
    vl1 <- viral_load_func(pars_tmp, ages)
    vl2 <- viral_load_func(pars_tmp_variant, ages)
    
    
    ## Modal Ct values
    modal_ct1 <- tibble(ct=vl1,age=ages,virus="Original variant")
    modal_ct2 <- tibble(ct=vl2,age=ages,virus="New variant")
    modal_ct_all[[i]] <- bind_rows(modal_ct1,modal_ct2) %>% mutate(virus=factor(virus,levels=variant_levels)) %>% mutate(samp=samp)
    
    ## Simulate Ct values
    sim_cts1 <- virosolver::simulate_viral_loads_example(ages,pars_tmp,N) %>% mutate(virus="Original variant")
    sim_cts2 <- virosolver::simulate_viral_loads_example(ages,pars_tmp_variant,N)%>% mutate(virus="New variant")
    sim_cts_all[[i]] <- bind_rows(sim_cts1, sim_cts2) %>% mutate(virus=factor(virus,levels=variant_levels)) %>% mutate(samp=samp)
  }
  
  modal_ct_all <- do.call("bind_rows",modal_ct_all)
  sim_cts_all <- do.call("bind_rows",sim_cts_all)
  
  modal_ct_summaries <- modal_ct_all %>% group_by(age, virus) %>% summarize(mean_ct=mean(ct),lower_50=quantile(ct,0.25),upper_50=quantile(ct,0.75),
                                                                            lower_95=quantile(ct,0.025),upper_95=quantile(ct,0.975))
  
  p <- ggplot(sim_cts_all %>% filter(ct < 40)) + 
    #geom_rect(ymin=42,ymax=vl_pars["intercept"],xmin=-1,xmax=max(ages)+1,fill="grey70",alpha=0.25) +
    #geom_jitter(aes(x=age,y=ct,col=virus),alpha=0.25,height=0,width=0.25,size=0.2) + 
    geom_ribbon(data=modal_ct_summaries,aes(x=age,ymin=lower_50,ymax=upper_50,fill=virus),alpha=0.5) +
    geom_ribbon(data=modal_ct_summaries,aes(x=age,ymin=lower_95,ymax=upper_95,fill=virus),alpha=0.1) +
    geom_line(data=modal_ct_summaries,aes(x=age,y=mean_ct,col=virus)) +
    coord_flip()+
    coord_cartesian(ylim=c(40,5),xlim=c(0,max(ages))) +
    scale_x_continuous(breaks=seq(0, max(ages),by=5)) +
    variant_color_scale + variant_fill_scale +
    xlab("Time since infection") + ylab("Ct value") +
    theme_overall + theme_nice_axes + theme(legend.position=c(0.7,0.8))
  
  if(!is.null(true_pars1)){
    p <- p +
      geom_line(data=true_modal_cts,aes(x=age,y=ct,col=virus),linetype="dashed")
  }
  
  return(p)
}


plot_virosolver_comparisons <- function(chain, true_vals, real_scales, real_v1_gr, real_v2_gr){
  chain$beta_diff <- chain$beta_alt - chain$beta
  chain_grs <- chain %>% dplyr::select(sampno, beta, beta_alt,beta_diff) %>% pivot_longer(-sampno)
  gr_key <- c("beta"="Original variant","beta_alt"="New variant", "beta_diff"="Difference")
  chain_grs$name <- gr_key[chain_grs$name]
  chain_grs$name <- factor(chain_grs$name,levels=c("Original variant","New variant","Difference"))
 
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
  
  p2 <- ggplot(chain_scales) + 
    geom_hline(yintercept=1,linetype="dashed")+
    geom_violin(aes(x=name,y=value),alpha=0.5,draw_quantiles=c(0.025,0.5,0.975),fill="red") + 
    geom_point(data=real_scales,aes(x=name,y=value)) +
    #scale_y_continuous(limits=c(0,7)) +
    scale_y_log10() +
    theme_overall +
    ylab("Relative value of new variant") +
    xlab("")
  
  
  v1_preds <- virosolver::plot_prob_infection(chain, 1000,exponential_growth_model, 0:35,true_prob_infection = real_v1_gr)
  v1_preds_quants <- v1_preds$predictions %>% group_by(t) %>% summarize(lower95=quantile(prob_infection,0.025),lower50=quantile(prob_infection,0.25),
                                                                        median=median(prob_infection),
                                                                        upper50=quantile(prob_infection,0.75),upper95=quantile(prob_infection,0.975)) 
  p_v1_preds <- ggplot(v1_preds_quants) + 
    geom_ribbon(aes(x=t,ymin=lower95,ymax=upper95),fill="blue",alpha=0.1) +
    geom_ribbon(aes(x=t,ymin=lower50,ymax=upper50),fill="blue",alpha=0.5) +
    geom_line(aes(x=t,y=median),col="blue") +
    geom_line(data=real_v1_gr,aes(x=t,y=prob_infection),linetype="dashed",col="black",size=0.75) +
    scale_y_continuous(limits=c(0,0.25)) +
    theme_overall +
    ylab("Original variant growth rate") +
    theme_no_x_axis
  
  chain1 <- chain
  chain1$beta <- chain1$beta_alt
  v2_preds <- virosolver::plot_prob_infection(chain1, 1000,exponential_growth_model, 0:35,true_prob_infection = real_v2_gr)
  v2_preds_quants <- v2_preds$predictions %>% group_by(t) %>% summarize(lower95=quantile(prob_infection,0.025),lower50=quantile(prob_infection,0.25),
                                                                        median=median(prob_infection),
                                                                        upper50=quantile(prob_infection,0.75),upper95=quantile(prob_infection,0.975)) 
  p_v2_preds <- ggplot(v2_preds_quants) + 
    geom_ribbon(aes(x=t,ymin=lower95,ymax=upper95),fill="red",alpha=0.1) +
    geom_ribbon(aes(x=t,ymin=lower50,ymax=upper50),fill="red",alpha=0.5) +
    geom_line(aes(x=t,y=median),col="red") +
    geom_line(data=real_v2_gr,aes(x=t,y=prob_infection),linetype="dashed",col="black",size=0.75) +
    scale_y_continuous(limits=c(0,0.25)) +
    theme_overall+
    ylab("New variant growth rate") +
    xlab("Time (35 days prior to sample)")
  
  v1_preds_quants$virus <- "Original variant"
  v2_preds_quants$virus <- "New variant"
  
  v_preds_comb <- bind_rows(v1_preds_quants,v2_preds_quants)
  
  real_grs_comb <- bind_rows(real_v1_gr %>% mutate(virus="Original variant"),
                             real_v2_gr %>% mutate(virus="New variant"))
  
  
  p_comb_preds <- ggplot(v_preds_comb) + 
    geom_ribbon(aes(x=t,ymin=lower95,ymax=upper95,fill=virus),alpha=0.1) +
    geom_ribbon(aes(x=t,ymin=lower50,ymax=upper50,fill=virus),alpha=0.5) +
    geom_line(aes(x=t,y=median,col=virus,linetype="Mean estimate")) +
    variant_color_scale_min + variant_fill_scale_min +
    geom_line(data=real_grs_comb,aes(x=t,y=prob_infection,col=virus,linetype="True value"),size=0.75) +
    scale_linetype_manual(name="",values=c("True value"="dashed","Mean estimate"="solid")) +
    #geom_line(data=real_v2_gr,aes(x=t,y=prob_infection),linetype="dashed",col="black",size=0.75) +
    theme_overall+
    theme_nice_axes +
    theme(legend.position=c(0.3,0.7),legend.text = element_text(size=6),legend.title = element_blank(),
          legend.spacing.y=unit(0,"cm"),legend.direction = "horizontal") +
    scale_x_continuous(breaks=seq(0,35,by=5),labels=rev(0-seq(0,35,by=5))) +
    ylab("Relative probability of infection") +
    xlab("Days (relative to sample date)")
  
  return(p_comb_preds)
}



plot_time_since_infection <- function(obs_times, vl_pars1, vl_pars2, prob_infection1, prob_infection2){
  variants <- c("Original variant","New variant")
  
  comb_dat <- NULL
  index <- 1
  for(j in seq_along(variants)){
    if(j == 1){
      pars <- vl_pars1
      prob_infection <- prob_infection1
    } else {
      pars <- vl_pars2
      prob_infection <- prob_infection2
    }
    max_age <- pars["max_incu_period"] + pars["max_sampling_delay"]
   
    ## Time at which standard deviation is reduced
    t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
    sd_mod <- rep(pars["sd_mod"], max_age) ## Up until t_switch, full standard deviation
    
    ## Prior to t_switch, sd=1
    unmod_vec <- 1:min(t_switch,max_age)
    sd_mod[unmod_vec] <- 1
    
    ## For the next sd_mod_wane days, variance about modal Ct trajectory decrease linearly
    decrease_vec <- (t_switch+1):(t_switch+pars["sd_mod_wane"])
    sd_mod[decrease_vec] <- 1 - ((1-pars["sd_mod"])/pars["sd_mod_wane"])*seq_len(pars["sd_mod_wane"])
    
    for(obs_time in obs_times){
      ## Restrict ages to times occurring before 
      ## time of sample collection (single cross section)
      
      ## Returns the full probability density distribution to simulate from
      densities <- pred_age_since_inf_symptomatic(pars["max_incu_period"],pars["max_sampling_delay"],obs_time,
                                                  pars, prob_infection, sd_mod)
      
      comb_dat[[index]] <- tibble(age=0:(max_age),density=densities, t=obs_time,virus=variants[j])
      index <- index + 1
    }
  }
  comb_dat <- do.call("bind_rows",comb_dat)

  prob_dist_1 <- comb_dat %>% filter(virus == "Original variant") %>% pull(density)
  prob_dist_2 <- comb_dat %>% filter(virus == "New variant") %>% pull(density)
  
  mean1 <- mean(sample(0:max_age,1000000,prob=prob_dist_1,replace=TRUE))
  mean2 <- mean(sample(0:max_age,1000000,prob=prob_dist_2,replace=TRUE))
  
  p <- ggplot(comb_dat) + 
    geom_ribbon(aes(x=age,ymin=0,ymax=density,fill=virus,group=virus),alpha=0.25,col="black") +
    geom_vline(data=tibble(virus=c("Original variant","New variant"), mean_delay=c(mean1,mean2)),aes(xintercept=mean_delay,col=virus),linetype="dashed") +
    variant_fill_scale_min + variant_color_scale_min +
    coord_cartesian(xlim=c(0,25)) +
    scale_x_continuous(breaks=seq(0,30,by=5)) +
    scale_y_continuous(expand=c(0,0),breaks=seq(0,1,by=0.025)) +
    ylab("Density") +
    xlab("Time since infection\n(incubation period + sampling delay)") +
    theme_overall + theme_nice_axes +
    theme(legend.position=c(0.8,0.8))
  p
}
