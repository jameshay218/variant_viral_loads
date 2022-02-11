## Calculate distribution skew based on PDF
calc_skew <- function(values,weights) {
  weights_std <- weights/sum(weights, na.rm=TRUE)
  xbar <- sum(values*weights_std, na.rm=TRUE)
  xi_xbar <- values - xbar
  return((sum(weights_std*xi_xbar^3))/((sum(weights_std*xi_xbar^2))^(3/2)))
}

## Calculate mode of distribution data
calc_mode <- function(values){
  if(length(values) < 2) return(NaN)
  dens <- density(values)
  return(dens$x[which.max(dens$y)])
}

## Calculate summary statistics from individual-level Ct simulations
calculate_ct_dist_indiv_summaries <- function(cts_indiv){
  cts_indiv %>% group_by(t, virus) %>%
    filter(ct < 40) %>%
    summarize(median=median(ct),mean=mean(ct),mode=calc_mode(ct),skewness=moments::skewness(ct),N=n())
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
    xlab("Time of epidemic") +
    theme_overall
  
  p_ages_relative <- ggplot(age_res) + 
    geom_tile(aes(x=t,y=infection_age,fill=scaled_proportion)) + 
    geom_line(data=age_res_mean,aes(x=t,y=mean_age),col="white",size=0.75) +
    scale_fill_viridis_c(name="Proportion of detectable infections") + 
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0),limits=c(0,max(ages))) +
    ylab("Days since infection") +
    xlab("Time of epidemic") +
    theme_overall
  
  p_prev <- ggplot(age_res %>% group_by(t) %>% mutate(prev=sum(proportion,na.rm=TRUE))) + 
    geom_line(aes(x=t,y=prev)) + 
    ylab("Prevalence") +
    xlab("Time of epidemic") +
    theme_overall
  
  return(list(age_distributions=age_res,p_ages=p_ages,p_ages_relative=p_ages_relative,p_prev=p_prev, age_mean=age_res_mean))  
}

## For a given incidence curve (per capita) and curve describing viral kinetics on the Ct scale
## Calculates the expected distribution of detectable Ct values and summary statistics for each day
calculate_ct_distribution <- function(vl_pars, ages, incidence, obs_times,symptom_surveillance=FALSE){
  ## Use virosolver to get the expected Ct distribution
  ct_distribution <- pred_dist_wrapper(seq(0,40,by=0.01),obs_times = obs_times,ages,vl_pars,incidence, symptom_surveillance)

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

simulate_cross_section <- function(pars, ages, incidence, obs_time,N=1000, use_pos=TRUE, symptom_surveillance=FALSE){
    ## Use virosolver to get the expected Ct distribution
    ct_distribution <- pred_dist_wrapper(seq(0,40,by=0.01),obs_times = obs_time,ages,pars,incidence,symptom_surveillance)
    ## Re-scale 
    if(use_pos){
      ct_distribution <- ct_distribution %>% filter(ct < 40)
    }
    
    ct_scaled_dist <- ct_distribution %>% group_by(t) %>% mutate(density_scaled=density/sum(density)) %>%  mutate(cumu_density=cumsum(density_scaled)) 
    cts <- sample(ct_scaled_dist$ct,size=N,replace=TRUE,prob=ct_scaled_dist$density_scaled)
    return(cts)
}


simulate_m_cross_sections <- function(pars, ages, incidence, obs_time,N=1000, m=10, use_pos=TRUE, symptom_surveillance=FALSE,align_gr=FALSE,grs=NULL,gr_test=NULL,virus_tmp=NULL){
  ## Use virosolver to get the expected Ct distribution
  if(is.null(dim(incidence))) {
    ct_distribution <- pred_dist_wrapper(seq(0,40,by=0.01),obs_times = obs_time,ages,pars,incidence,symptom_surveillance)
    ## Re-scale 
    if(use_pos){
      ct_distribution <- ct_distribution %>% filter(ct < 40)
    }
    
    ct_scaled_dist <- ct_distribution %>% group_by(t) %>% mutate(density_scaled=density/sum(density)) %>%  mutate(cumu_density=cumsum(density_scaled)) 
  }
  
  ## Sample from distribution m times and combine
  cts <- list()
  for(j in 1:m){
    if(!is.null(dim(incidence))) {
      run_no <- sample(unique(incidence$run),1)
      if(align_gr){
        samp_times <- grs %>% 
          filter(inc > 0.0001) %>%
          group_by(virus) %>% 
          mutate(gr_diff=abs(gr - gr_test)) %>% 
          filter(gr_diff==min(gr_diff,na.rm=TRUE),
                 virus==virus_tmp)
        #print(samp_times)
        obs_time <- samp_times %>% pull(t)
      } 
      
      incidence_tmp <- incidence %>% dplyr::filter(run == run_no) %>% pull(inc)
      
      ct_distribution <- pred_dist_wrapper(seq(0,40,by=0.01),obs_times = obs_time,ages,pars,incidence_tmp,symptom_surveillance)
      ## Re-scale 
      if(use_pos){
        ct_distribution <- ct_distribution %>% filter(ct < 40)
      }
      
      ct_scaled_dist <- ct_distribution %>% group_by(t) %>% mutate(density_scaled=density/sum(density)) %>%  mutate(cumu_density=cumsum(density_scaled)) 
    } 
    
    
    cts[[j]] <- tibble(ct=sample(ct_scaled_dist$ct,size=N,replace=TRUE,prob=ct_scaled_dist$density_scaled),sampno=j,obs_time=obs_time)
  }
  cts <- do.call("bind_rows",cts)
  return(cts)
}



## Simulate individual-level kinetics from the pooling paper model, passing in the MCMC chain
simulate_individual_level_data <- function(N, incidence, times, chain, parTab, max_vl=11,ct_intercept=40,vl_lod=2,obs_time=NULL){
  
  samp_times <- times
  if(!is.null(obs_time)) {
    samp_times  <- obs_times
  }
  
  ## Simulate infection times
  infection_times <- simulate_infection_times(N, sum(incidence),incidence)
  infection_times_dat <- tibble(i=seq_along(infection_times), inf_time=infection_times)
  
  ## Simulate viral loads for the sample population
  ## <0.2% of simulated viral loads are 11 or higher
  simulated_data <- simulate_viral_loads_hinge(infection_times, samp_times, chain, parTab,save_during=FALSE,
                                               add_noise=TRUE,max_vl=max_vl,simno=NA)
  
  
  viral_loads <- simulated_data$viral_loads
  viral_loads_melted <- reshape2::melt(viral_loads)
  colnames(viral_loads_melted) <- c("i","t","viral_load")
  viral_loads_melted <- as_tibble(viral_loads_melted) %>% left_join(infection_times_dat)
  
  obs_dat <- simulated_data$obs
  obs_vl_melted <- reshape2::melt(obs_dat)
  colnames(obs_vl_melted) <- c("i","t","obs")
  obs_vl_melted <- as_tibble(obs_vl_melted) %>% left_join(infection_times_dat)
  
  sim_pars <- as.data.frame(simulated_data$pars)
  colnames(sim_pars)[5] <- "incu_period"
  sim_pars$i <- 1:nrow(sim_pars)
  sim_pars <- sim_pars %>% as_tibble() %>% dplyr::select(i, incu_period)
  
  obs_vl_melted <- obs_vl_melted %>% left_join(sim_pars) %>% left_join(viral_loads_melted) %>% left_join(infection_times_dat) 
  obs_vl_melted <- obs_vl_melted %>% mutate(days_since_infection = t - inf_time) %>% arrange(i, days_since_infection) %>% 
    mutate(ct=ct_intercept - log2(10)*(obs-vl_lod)) %>%
    mutate(ct=ifelse(ct < 0,0, ct),obs=ifelse(obs > max_vl, max_vl, obs))
  
  return(obs_vl_melted)
}



## Simulate individual-level kinetics from the pooling paper model, passing in the MCMC chain
simulate_individual_level_data_symptomatic <- function(N, incidence, times, chain, parTab, max_vl=11,ct_intercept=40,vl_lod=2,confirm_delays){
  
  ## Simulate infection times
  infection_times <- simulate_infection_times(N, sum(incidence),incidence)
  infection_times_dat <- tibble(i=seq_along(infection_times), inf_time=infection_times)
  
  ## Simulate viral loads for the sample population
  ## <0.2% of simulated viral loads are 11 or higher
  simulated_data <- simulate_viral_loads_hinge(infection_times, samp_times, chain, parTab,save_during=FALSE,
                                               add_noise=TRUE,max_vl=max_vl,simno=NA,
                                               symptom_surveillance = TRUE, confirm_delays=confirm_delays)
  
  
  viral_loads <- simulated_data$viral_loads
  viral_loads_melted <- reshape2::melt(viral_loads)
  colnames(viral_loads_melted) <- c("i","t","viral_load")
  viral_loads_melted <- as_tibble(viral_loads_melted) %>% left_join(infection_times_dat)
  
  obs_dat <- simulated_data$obs
  obs_vl_melted <- reshape2::melt(obs_dat)
  colnames(obs_vl_melted) <- c("i","t","obs")
  obs_vl_melted <- as_tibble(obs_vl_melted) %>% left_join(infection_times_dat)
  
  sim_pars <- as.data.frame(simulated_data$pars)
  colnames(sim_pars)[5] <- "incu_period"
  sim_pars$i <- 1:nrow(sim_pars)
  sim_pars <- sim_pars %>% as_tibble() %>% dplyr::select(i, incu_period)
  
  obs_vl_melted <- obs_vl_melted %>% left_join(sim_pars) %>% left_join(viral_loads_melted) %>% left_join(infection_times_dat) 
  obs_vl_melted <- obs_vl_melted %>% mutate(days_since_infection = t - inf_time) %>% arrange(i, days_since_infection) %>% 
    mutate(ct=ct_intercept - log2(10)*(obs-vl_lod))
  
  obs_vl_melted$confirm_delay <- confirm_delays
  
  obs_vl_melted <- obs_vl_melted %>% mutate(obs_time = round(inf_time + incu_period+confirm_delay))
  
  return(obs_vl_melted)
}

resample_ct_dist <- function(obs, samp_size=100, N=1,with_replacement=FALSE,bootstrap_cts=TRUE,cts=seq(0,40,by=0.1),bde_bandwidth=0.03){
  n_samp <- samp_size
  if(with_replacement==FALSE){
    n_samp <- min(length(obs), samp_size)
  }
  ## Each column is one trial, each row is one sample per trial
  if(bootstrap_cts){
    samp_cts <- sapply(1:N, function(x) sample(obs, size=n_samp,replace=with_replacement))
  } else {
    dens_estimate1 <- bde::bde(obs,estimator="betakernel",upper.limit = 40,lower.limit=0,b=bde_bandwidth)
    
    #dens_estimate2 <- bde::bde(obs,estimator="vitale",upper.limit = 39.9,lower.limit=0,b=bde_bandwidth)
    #dens_estimate3 <- bde::bde(obs,estimator="boundarykernel",upper.limit = 39.9,lower.limit=0,b=bde_bandwidth)
    #dens_estimate4 <- bde::bde(obs,estimator="kakizawa",upper.limit = 39.9,lower.limit=0,b=bde_bandwidth)
    
    dens1 <- bde::density(dens_estimate1,cts)
    
    #dens2 <- density(dens_estimate2,cts)
    #dens3 <- density(dens_estimate3,cts)
    #dens4 <- density(dens_estimate4,cts)
    
    #plot(cts,dens1,type='l')
    #lines(cts, dens2,col="red")
    #lines(cts, dens3,col="blue")
    #lines(cts, dens4,col="green")
    
    ## Each column is one trial, each row is one sample per trial
    samp_cts <- sapply(1:N, function(x) sample(cts, size=n_samp,replace=with_replacement,prob = dens1))
  }
  samp_cts <- reshape2::melt(samp_cts)
  colnames(samp_cts) <- c("i","trial","ct")

  return(samp_cts)
}


create_posterior_func_compare <- function(parTab,
                                          data,
                                          PRIOR_FUNC=NULL,
                                          INCIDENCE_FUNC=NULL,
                                          solve_ver="likelihood",
                                          solve_likelihood=TRUE,
                                          use_pos=FALSE,
                                          symptom_surveillance=FALSE,
                                          ...) {
  par_names <- parTab$names
  pars <- parTab$values
  names(pars) <- par_names
  times <- 0:max(data$t)
  ages <- 1:max(data$t)
  obs_times <- unique(data$t)
  variants <- unique(data$variant)
  
  ## Pull out undetectable Ct values and count how many per observation time
  ## We only need get the likelihood of an undetectable Ct value once per time point,
  ## and then just have this contribute N times where N is the number of undetectable Cts
  ## at that time point
  data_use <- data
  undetectable_counts <- NULL
  if("intercept" %in% par_names){
    undetectable_tally <- data_use %>%
      filter(ct >= pars["intercept"]) %>%
      group_by(t) %>%
      tally()
    no_undetectable_times <- setdiff(obs_times, unique(undetectable_tally$t))
    no_undetectable_tally <- tibble(t=no_undetectable_times,n=0)
    undetectable_tally <- bind_rows(undetectable_tally, no_undetectable_tally) %>% arrange(t)
    undetectable_counts <- undetectable_tally$n
    data_use <- data_use %>% filter(ct < pars["intercept"])
  }
  
  ## Pull out data into a list for quicker indexing later on
  data_list <- NULL
  for(variant1 in variants){
    data_list[[variant1]] <- list()
    for(i in seq_along(obs_times)){
      data_list[[variant1]][[i]] <- data_use %>% filter(t == obs_times[i],variant==variant1) %>% pull(ct)
    }
  }
  
  
  f <- function(pars){
    ## For each variant
    lik <- 0
    preds <- NULL
    pars_list <- NULL
    prob_infection_list <- NULL
    
    for(variant1 in variants){
      pars_tmp <- pars
      names(pars_tmp) <- par_names
      if(variant1 != variants[1]){
        ## Growth rate pars
        #pars_tmp["beta"] <- pars_tmp["beta"]*pars_tmp["beta_scale"]
        pars_tmp["beta"] <- pars_tmp["beta_alt"]
        pars_tmp["overall_prob"] <- pars_tmp["overall_prob_scale"]
        
        if("prob" %in% par_names){
          pars_tmp[which(par_names == "prob")] <- pars_tmp[which(par_names == "prob_alt")]
        }
        
        ## Viral load control points
        pars_tmp["viral_peak"] <- pars_tmp["viral_peak"]*pars_tmp["viral_peak_scale"]
        pars_tmp["level_switch"] <- pars_tmp["level_switch"]*pars_tmp["level_switch_scale"]
        
        ## Timing pars
        pars_tmp["desired_mode"] <- pars_tmp["desired_mode"]*pars_tmp["desired_mode_scale"]
        pars_tmp["t_switch"] <- pars_tmp["t_switch"]*pars_tmp["t_switch_scale"]
        
        ## Detect loss
        pars_tmp["prob_detect"] <- pars_tmp["prob_detect"]*pars_tmp["prob_detect_scale"]
        
        ## Measurement parameters
        pars_tmp["obs_sd"] <- pars_tmp["obs_sd"]*pars_tmp["obs_sd_scale"]
        pars_tmp["sd_mod"] <- pars_tmp["sd_mod"]*pars_tmp["sd_mod_scale"]
        pars_tmp["sd_mod_wane"] <- pars_tmp["sd_mod_wane"]*pars_tmp["sd_mod_wane_scale"]
        
      }
      prob_infection_tmp <- INCIDENCE_FUNC(pars_tmp, times)
      
      pars_list[[variant1]] <- pars_tmp
      prob_infection_list[[variant1]] <- prob_infection_tmp
      
      if(solve_ver == "likelihood"){
        if(solve_likelihood){
          lik <- lik + sum(likelihood_cpp_wrapper(data_list[[variant1]], ages, obs_times,pars_tmp, prob_infection_tmp,
                                                  use_pos,0))
        }
        
        if(!is.null(PRIOR_FUNC)){
          prior <- PRIOR_FUNC(pars_tmp, ...)
          lik <- lik+prior
        }
      } else {
        preds <- bind_rows(preds, pred_dist_wrapper(seq(0,40,by=1),obs_times,ages,pars,prob_infection_tmp,symptom_surveillance) %>% mutate(variant=variant1))
      }
    }
    if(solve_ver == "likelihood"){
      ## Now solve undetectable likelihood
      if(!use_pos){
        lik <- lik + likelihood_undetectable_2variants(undetectable_counts, variants,obs_times, ages,
                                                       pars_list, prob_infection_list)
      }
    return(lik) 
    } else {
      return(preds)
    }
    f
  }
}

likelihood_undetectable_2variants <- function(undetectable_counts, variants,times, ages, pars_list, 
                                              prob_infection_list){
  liks <- 0
  for(i in seq_along(times)){
    n_obs <- undetectable_counts[i]
    p_undetectables <- numeric(length(variants))
    
    for(j in seq_along(variants)){
      variant1 <- variants[[j]]
      pars_tmp <- pars_list[[j]]
      prob_infection_tmp <- prob_infection_list[[j]]

      ## Because we have a different standard deviation for different times
      ## Time at which standard deviation is reduced
      t_switch <-  pars_tmp["t_switch"] + pars_tmp["desired_mode"] + pars_tmp["tshift"]
      sd_mod <- rep(pars_tmp["sd_mod"], max(ages))
      
      ## Prior to t_switch, full standard deviation
      ## Make sure we don't go past maximum of ages vector
      unmod_vec <- 1:min(t_switch,max(ages)+1)
      sd_mod[unmod_vec] <- 1
      
      ## For the next sd_mod_wane days, decrease linearly
      decrease_vec <- (t_switch+1):(t_switch+pars_tmp["sd_mod_wane"])
      ## The rest are at sd_mod
      sd_mod[decrease_vec] <- 1 - ((1-pars_tmp["sd_mod"])/pars_tmp["sd_mod_wane"])*seq_len(pars_tmp["sd_mod_wane"])
      
      ## Only solve back until the earliest possible infection time
      ages1 <- ages[(times[i] - ages) > 0]
      sd_mod1 <- sd_mod[(times[i] - ages) > 0]
      
      ## Find the log probability of being undetectable from infection with this strain,
      ## exponential and subtract from 1 to give the probability of being detectable and
      ## infected with this variant
      p_undetectables[j] <-  1-exp(likelihood_cpp(pars_tmp["intercept"], times[i], ages1,
                                                  pars_tmp, prob_infection_tmp,
                                          sd_mod))
    }
    ## Find the probability of not being detectable and infected with any variant
    p_undetectable <- 1 - sum(p_undetectables)
    
    ## Log and multiply to get log probability
    liks <- liks + n_obs*log(p_undetectable)
  }
  liks
}


