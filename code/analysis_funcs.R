## Calculate distribution skew based on PDF
calc_skew <- function(values,weights) {
  weights_std <- weights/sum(weights, na.rm=TRUE)
  xbar <- sum(values*weights_std, na.rm=TRUE)
  xi_xbar <- values - xbar
  return((sum(weights_std*xi_xbar^3))/((sum(weights_std*xi_xbar^2))^(3/2)))
}

## For a given incidence curve (per capita) and curve describing the proportion of individuals remaining 
## PCR detectable on each day post infection, calculates the proportion of individuals of each infection ages
## in days over the course of the epidemic
calculate_infection_age_distribution <- function(incidence, detectable_props, ages){
  ## Find distribution of ages since infection among PCR positive individuals as convolution of incidence curve and proportion detectable
  lastday <- length(detectable_props)
  n_times <- length(incidence)
  age_res <- matrix(nrow=n_times,ncol=lastday)
  for (i in 2:n_times) {
    past_inc <- incidence[(i-1):(max(i-lastday,1))]
    days <- seq_along(past_inc)
    age_res[i,days] <- past_inc*detectable_props[days]
  }
  
  ## Get mean age of infection
  age_res_overall <- age_res/apply(age_res, 1, sum, na.rm=TRUE)
  age_mean <- tibble(t=seq_along(incidence),mean_age=apply(age_res_overall, 1, function(res) sum(res*(1:lastday), na.rm=TRUE))) %>% filter(t > lastday)
  
  age_res <- as_tibble(age_res)
  colnames(age_res) <- ages
  age_res <- age_res %>% mutate(t=seq_along(incidence))
  age_res <- age_res %>% pivot_longer(-t,names_to="infection_age",names_transform=list(infection_age=as.integer),values_to="proportion") %>%
    group_by(t) %>% mutate(prev=sum(proportion), scaled_proportion=proportion/prev) %>% ungroup()
  age_res_mean <- age_res %>% group_by(t) %>% mutate(proportion_weights=scaled_proportion*infection_age) %>% summarize(mean_age=sum(proportion_weights))
  
  p_ages <- ggplot(age_res) + 
    geom_tile(aes(x=t,y=infection_age,fill=proportion)) + 
    scale_fill_viridis_c(name="Proportion of all individuals") + 
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0),limits=c(0,max(ages))) +
    ylab("Days since infection") +
    xlab("Time of epidemice") +
    theme_overall
  
  p_ages_relative <- ggplot(age_res) + 
    geom_tile(aes(x=t,y=infection_age,fill=scaled_proportion)) + 
    geom_line(data=age_res_mean,aes(x=t,y=mean_age),col="white",size=0.75) +
    scale_fill_viridis_c(name="Proportion of detectable infections") + 
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0),limits=c(0,max(ages))) +
    ylab("Days since infection") +
    xlab("Time of epidemice") +
    theme_overall
  
  p_prev <- ggplot(age_res %>% group_by(t) %>% mutate(prev=sum(proportion,na.rm=TRUE))) + 
    geom_line(aes(x=t,y=prev)) + 
    ylab("Prevalence") +
    xlab("Time of epidemice") +
    theme_overall
  
  return(list(age_distributions=age_res,p_ages=p_ages,p_ages_relative=p_ages_relative,p_prev=p_prev, age_mean=age_res_mean))  
}

## For a given incidence curve (per capita) and curve describing viral kinetics on the Ct scale
## Calculates the expected distribution of detectable Ct values and summary statistics for each day
calculate_ct_distribution <- function(vl_pars, ages, incidence, obs_times){
  ## Use virosolver to get the expected Ct distribution
  ct_distribution <- pred_dist_wrapper(seq(0,40,by=0.01),obs_times = obs_times,ages,vl_pars,incidence)
  
  ## Re-scale 
  ct_scaled_dist <- ct_distribution %>% filter(ct < 40) %>% group_by(t) %>% mutate(density_scaled=density/sum(density)) %>%  mutate(cumu_density=cumsum(density_scaled)) 
  
  ## Calculate skew of detectable Cts over time
  ct_skew <- ct_scaled_dist %>% filter(ct < 40) %>% group_by(t) %>% summarize(skewness=calc_skew(ct,density_scaled))
  ## Calculate mean Ct based on weighted distribution
  ct_mean <- ct_scaled_dist %>% summarize(mean=sum(ct*density_scaled))
  ## Use cumulative density to get median and central 50% quantiles
  ct_med <- ct_scaled_dist %>% filter(cumu_density >= 0.5) %>% filter(ct == min(ct))%>% dplyr::select(ct, t) %>% rename(median=ct)
  ct_lower <- ct_scaled_dist %>% filter(cumu_density >= 0.25) %>% filter(ct == min(ct)) %>% dplyr::select(ct, t) %>% rename(low25=ct)
  ct_upper <- ct_scaled_dist %>% filter(cumu_density >= 0.75) %>% filter(ct == min(ct))%>% dplyr::select(ct, t) %>% rename(upp75=ct)
  ct_summary <- left_join(ct_med,ct_lower) %>% left_join(ct_upper) %>% left_join(ct_mean) %>% left_join(ct_skew)
  return(ct_summary)
}

simulate_cross_section <- function(pars, ages, incidence, obs_time,N=1000){
    ## Use virosolver to get the expected Ct distribution
    ct_distribution <- pred_dist_wrapper(seq(0,40,by=0.01),obs_times = obs_time,ages,pars,incidence)
    ## Re-scale 
    ct_scaled_dist <- ct_distribution %>% filter(ct < 40) %>% group_by(t) %>% mutate(density_scaled=density/sum(density)) %>%  mutate(cumu_density=cumsum(density_scaled)) 
    cts <- sample(ct_scaled_dist$ct,size=N,replace=TRUE,prob=ct_scaled_dist$density_scaled)
    return(cts)
}
