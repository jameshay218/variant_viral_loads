theme_no_x_axis <- theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank(),axis.line.x=element_blank())
theme_nice_axes <- theme(axis.line.x=element_line(size=0.5,color="black"),axis.line.y=element_line(size=0.5,color="black"),panel.border = element_blank())
theme_overall <- theme_bw() + theme(axis.text=element_text(size=7),axis.title = element_text(size=8),plot.tag = element_text(size=10,face="bold"),
                                    legend.text=element_text(size=7),legend.title=element_text(size=8))

variant_levels <- c("Overall","Original variant","New variant","New variant, same kinetics","New variant, different kinetics")

color_key <- c("Overall"="black",
               "Original variant"="blue",
               "New variant"="red",
               "New variant, same kinetics"="red",
               "New variant, different kinetics"="orange")

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

variant_color_scale <- scale_color_manual(name="Variant",values=color_key) 
variant_fill_scale <- scale_fill_manual(name="Variant",values=color_key) 
variant_linetype_scale <- scale_linetype_manual(name="Variant",values=linetype_key)
variant_size_scale <- scale_size_manual(name="Variant",values=size_key)

plot_medians_and_skew <- function(combined_summaries){
  p_medians <- ggplot(combined_summaries) + 
    geom_line(aes(x=t,y=median,col=virus,linetype=virus,size=virus)) +
    variant_color_scale + variant_fill_scale + variant_linetype_scale + variant_size_scale +
    scale_y_continuous(trans="reverse") +
    scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
    ylab("Median Ct") +
    theme_overall +
    xlab("Time") +
    theme(legend.position="none") +theme_nice_axes# + theme_no_x_axis
  
  p_skews <- ggplot(combined_summaries) + 
    geom_line(aes(x=t,y=skewness,col=virus,linetype=virus,size=virus))+
    variant_color_scale + variant_fill_scale + variant_linetype_scale + variant_size_scale +
    scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
    ylab("Skewness of Cts") +
    theme_overall +
    xlab("Time") +
    theme(legend.position="none") + theme_nice_axes
  
  return(list(p_medians,p_skews))
}

plot_growth_rate_lineups <- function(ct_combined_summaries){
  p_medians <- ct_combined_summaries %>% pivot_longer(-c(t,virus,gr,inc)) %>% filter(name %in% c("median")) %>% filter(virus != "Overall") %>%
    ggplot() +
    geom_line(aes(x=gr,y=value,col=virus,linetype=virus,size=virus)) +
    ylab("Median Ct") +
    xlab("Growth rate") +
    scale_x_continuous(trans="reverse") +
    scale_y_continuous(trans="reverse") +
    variant_color_scale + variant_linetype_scale + variant_size_scale +
    theme_overall + theme_nice_axes + theme(legend.position="none")
  
  p_skews <- ct_combined_summaries %>% pivot_longer(-c(t,virus,gr,inc)) %>% filter(name %in% c("skewness")) %>% filter(virus != "Overall") %>%
    ggplot() +
    geom_line(aes(x=gr,y=value,col=virus,linetype=virus,size=virus)) +
    scale_x_continuous(trans="reverse") +
    ylab("Skewness of Cts") +
    xlab("Growth rate") +
    variant_color_scale +
    theme_overall + variant_linetype_scale + variant_size_scale +
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
    geom_line(data=tibble(x=ages,y=modal_ct),aes(x=x,y=y),col="blue") +
    geom_line(data=tibble(x=ages,y=modal_ct),aes(x=x,y=y),col="red",linetype="dashed") +
    coord_flip()+
    coord_cartesian(ylim=c(40,5),xlim=c(0,max(ages))) +
    scale_x_continuous(breaks=seq(0, max(ages),by=5)) +
    xlab("Time since infection") + ylab("Ct value") +
    theme_overall + theme_nice_axes
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
    variant_color_scale +
    scale_linetype_manual(name="Variant",values=c("Original variant"="solid","New variant, same kinetics"="dashed", "New variant, different kinetics"="solid")) +
    xlab("Time since infection") + ylab("Ct value") +
    theme_overall + theme_nice_axes + theme(legend.position=c(0.7,0.8))
  return(p)
}


p_sim_ct_compare_naive <- function(vl_pars1,vl_pars2,virus1_inc,virus2_inc, ages,samp_time,N=100,dotsize=1) {
  cts_1 <- tibble(ct=simulate_cross_section(vl_pars1, ages, virus1_inc,obs_time=samp_time,N=N),virus="Original variant")
  cts_2 <- tibble(ct=simulate_cross_section(vl_pars1, ages, virus2_inc,obs_time=samp_time,N=N),virus="New variant, same kinetics")
  cts_2_alt <- tibble(ct=simulate_cross_section(vl_pars2, ages, virus2_inc,obs_time=samp_time,N=N),virus="New variant, different kinetics")
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



p_sim_ct_compare_growth <- function(vl_pars1,vl_pars2,virus1_inc,virus2_inc, ages,combined_summaries,growth_rate_samp,N=100,dotsize=1) {
  samp_times <- ct_combined_summaries %>% 
    group_by(virus) %>% 
    mutate(gr_diff=abs(gr - growth_rate_samp)) %>% 
    filter(gr_diff==min(gr_diff))
  
  samp_time1 <- samp_times %>% filter(virus == "Original variant") %>% pull(t)
  samp_time2 <- samp_times %>% filter(virus == "New variant, same kinetics") %>% pull(t)
  samp_time2_alt <- samp_times %>% filter(virus == "New variant, different kinetics") %>% pull(t)
  
  cts_1 <- tibble(ct=simulate_cross_section(vl_pars1, ages, virus1_inc,obs_time=samp_time1,N=N),virus="Original variant")
  cts_2 <- tibble(ct=simulate_cross_section(vl_pars1, ages, virus2_inc,obs_time=samp_time2,N=N),virus="New variant, same kinetics")
  cts_2_alt <- tibble(ct=simulate_cross_section(vl_pars2, ages, virus2_inc,obs_time=samp_time2_alt,N=N),virus="New variant, different kinetics")
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
