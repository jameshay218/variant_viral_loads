simulate_popn_cts_symptomatic <- function(virus1_inc, virus2_inc, 
                                          vl_pars1, vl_pars2,
                                          population_n, times,
                                          confirm_delay_par1_v1, confirm_delay_par2_v1,
                                          confirm_delay_par1_v2, confirm_delay_par2_v2,
                                          incu_period_par1_v1=1.621, incu_period_par2_v1=0.418,
                                          incu_period_par1_v2=1.621, incu_period_par2_v2=0.418){
  
  ###############################################
  ## STOCHASTIC SIMULATION OF CT VALUES
  ###############################################
  ## Simulate complete line list for individuals infection with the original or new variant
  v1_linelist <- virosolver::simulate_observations_wrapper(floor(virus1_inc*population_n),times=times,population_n=population_n,
                                                           incu_period_par1=incu_period_par1_v1, incu_period_par2=incu_period_par2_v1,
                                                           conf_delay_par1 = confirm_delay_par1_v1, conf_delay_par2 = confirm_delay_par2_v1) %>% 
    mutate(virus="Original variant") %>% filter(!is.na(infection_time)) %>% mutate(i=1:n())
  v2_linelist <- virosolver::simulate_observations_wrapper(floor(virus2_inc*population_n),times=times,population_n=population_n,
                                                           incu_period_par1=incu_period_par1_v2, incu_period_par2=incu_period_par2_v2,
                                                           conf_delay_par1 = confirm_delay_par1_v2, conf_delay_par2 = confirm_delay_par2_v2) %>% 
    mutate(virus="New variant") %>% filter(!is.na(infection_time))%>% mutate(i=1:n())
  v2_linelist$i <- v2_linelist$i + max(v1_linelist$i)
  
  ## Simulate situation where all individuals are observed at some point
  v1_observed_linelist <- simulate_reporting(v1_linelist %>% filter(is_infected == 1),frac_report=1,timevarying_prob=NULL,solve_times=times,symptomatic=TRUE)
  v2_observed_linelist <- simulate_reporting(v2_linelist %>% filter(is_infected == 1),frac_report=1,timevarying_prob=NULL,solve_times=times,symptomatic=TRUE)
  
  print("Simulating viral loads")
  ## Simulate Ct values from random cross sections
  cts_pop_sympt_v1 <- simulate_viral_loads_wrapper(v1_observed_linelist$sampled_individuals %>% filter(!is.na(infection_time)),kinetics_pars=vl_pars1)
  print("Virus 1 done")
  cts_pop_sympt_v2 <- simulate_viral_loads_wrapper(v2_observed_linelist$sampled_individuals %>% filter(!is.na(infection_time)),kinetics_pars=vl_pars1) %>% mutate(virus="New variant, same kinetics")
  print("Virus 2a done")
  cts_pop_sympt_v2_alt <- simulate_viral_loads_wrapper(v2_observed_linelist$sampled_individuals %>% filter(!is.na(infection_time)),kinetics_pars=vl_pars2) %>% mutate(virus="New variant, different kinetics")
  print("Virus 2b done")
  
  cts_pop_sympt_comb <- bind_rows(cts_pop_sympt_v1,cts_pop_sympt_v2,cts_pop_sympt_v2_alt) %>% mutate(virus=factor(virus,levels=variant_levels))
  
  cts_pop_summaries <- cts_pop_sympt_comb %>% filter(ct_obs < 40) %>% group_by(sampled_time,virus) %>% 
    summarize(mean_ct=mean(ct_obs),median_ct=median(ct_obs),
              mode_ct=calc_mode(ct_obs),skew_ct=moments::skewness(ct_obs),
              N=n())
  
  p_mean_sympt <- ggplot() + 
    geom_smooth(data=cts_pop_sympt_comb %>% filter(ct_obs < 40) %>% ungroup() %>% sample_n(min(n(),100000)),aes(x=sampled_time,y=ct_obs,col=virus,fill=virus),alpha=0.1) +
    variant_color_scale + variant_fill_scale +
    scale_y_continuous(trans="reverse") +
    scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
    theme_overall +
    theme(legend.position="none") +theme_nice_axes + theme_no_x_axis +
    xlab("Time") +
    ylab("Smoothed detectable Ct values") 
  
  return(list(
    ct_dat_sympt=cts_pop_sympt_comb, ct_dat_summary=cts_pop_summaries,
    p_mean_sympt=p_mean_sympt))
  
}