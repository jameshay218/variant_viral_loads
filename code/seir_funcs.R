## Takes a vector of model parameters, initial states and solve times and solves the two-strain SEIR model
run_2strain_seir_simulation <- function(pars, states, times, model_ver="lsoda",
                                        seir_filename="../code/odin_files/odinseirstrains.R",
                                        n_repeats=1,n_threads=4,seed=1){
  if(model_ver == "odin"){
    seir_sim <- simulate_odin_2strain_seir(seir_filename, max(times), pars, 
                                           n_threads=n_threads,
                                           seed=seed,n_repeats=n_repeats)
    
    virus1_inc_repeats_determ <- seir_sim$dat %>%
      dplyr::filter(strain %in% c("Original variant","Reinfection (new->original)"),compartment=="obs_lambda") %>% 
      group_by(time,run) %>% dplyr::summarize(inc=sum(y)) %>%
      dplyr::mutate(inc=inc) %>% #/pars["S_ini"]) %>% 
      dplyr::select(run,time,inc)
    
    virus2_inc_repeats_determ  <- seir_sim$dat %>%
      dplyr::filter(strain %in% c("New variant","Reinfection (original->new)"),compartment=="obs_lambda") %>% 
      group_by(time,run) %>% dplyr::summarize(inc=sum(y)) %>%
      dplyr::mutate(inc=inc) %>%#/pars["S_ini"]) %>% 
      dplyr::select(run,time,inc)
    
   virus1_inc_determ  <- seir_sim$dat %>%
      dplyr::filter(strain %in% c("Original variant","Reinfection (new->original)"),compartment=="obs_lambda") %>% 
      group_by(time,run) %>% dplyr::summarize(inc=sum(y)) %>%
      dplyr::mutate(inc=inc) %>% #/pars["S_ini"]) %>% 
      group_by(time) %>%
      dplyr::summarize(inc=mean(inc)) %>%
      pull(inc)
    
    virus2_inc_determ  <- seir_sim$dat %>%
      dplyr::filter(strain %in% c("New variant","Reinfection (original->new)"),compartment=="obs_lambda") %>% 
      group_by(time,run) %>% dplyr::summarize(inc=sum(y)) %>%
      dplyr::mutate(inc=inc) %>%#/pars["S_ini"]) %>% 
      group_by(time) %>%
      dplyr::summarize(inc=mean(inc)) %>%
      pull(inc)
    
    
    
    virus1_inc_repeats <- seir_sim$dat %>%
      dplyr::filter(strain %in% c("Original variant","Reinfection (new->original)"),compartment=="obs_inc") %>% 
      group_by(time,run) %>% dplyr::summarize(inc=sum(y)) %>%
      dplyr::mutate(inc=inc) %>% #/pars["S_ini"]) %>% 
      dplyr::select(run,time,inc)
    
    virus2_inc_repeats <- seir_sim$dat %>%
      dplyr::filter(strain %in% c("New variant","Reinfection (original->new)"),compartment=="obs_inc") %>% 
      group_by(time,run) %>% dplyr::summarize(inc=sum(y)) %>%
      dplyr::mutate(inc=inc) %>%#/pars["S_ini"]) %>% 
      dplyr::select(run,time,inc)
    
    virus1_inc <- seir_sim$dat %>%
      dplyr::filter(strain %in% c("Original variant","Reinfection (new->original)"),compartment=="obs_inc") %>% 
      group_by(time,run) %>% dplyr::summarize(inc=sum(y)) %>%
      dplyr::mutate(inc=inc) %>% #/pars["S_ini"]) %>% 
      group_by(time) %>%
      dplyr::summarize(inc=mean(inc)) %>%
      pull(inc)
    
    virus2_inc <- seir_sim$dat %>%
      dplyr::filter(strain %in% c("New variant","Reinfection (original->new)"),compartment=="obs_inc") %>% 
      group_by(time,run) %>% dplyr::summarize(inc=sum(y)) %>%
      dplyr::mutate(inc=inc) %>%#/pars["S_ini"]) %>% 
      group_by(time) %>%
      dplyr::summarize(inc=mean(inc)) %>%
      pull(inc)
    
    
    p_inc <- seir_sim$plot
    p_compartments <- NULL
    seir_sim <- seir_sim$dat
    
  } else {
    seir_sim <- as_tibble(as.data.frame(deSolve::lsoda(states,times,seir_model_2strains,pars)))
    ## Sense check compartments
    p_compartments <- seir_sim %>% 
      pivot_longer(-time) %>% 
      filter(!(name %in% c("inc1","inc2"))) %>% 
      ggplot() + geom_line(aes(x=time,y=value,col=name)) +
      theme_overall
    
    
    ## EXTRACT PER CAPITA INCIDENCE FOR EACH VIRUS
    virus1_inc <- c(0,diff(seir_sim$inc1))
    virus2_inc <- c(0,diff(seir_sim$inc2))
    
    virus1_inc_repeats <- NULL
    virus2_inc_repeats <- NULL
    
    virus1_inc_determ <- NULL
    virus2_inc_determ <- NULL
    virus1_inc_repeats_determ <- NULL
    virus2_inc_repeats_determ <- NULL
    
    seir_inc_dat1 <- tibble(t=times,inc=virus1_inc,virus="Original variant")
    seir_inc_dat2 <- tibble(t=times,inc=virus2_inc,virus="New variant")
    seir_inc_dat_comb <- tibble(t=times,inc=virus1_inc+virus2_inc,virus="Overall")
    seir_inc_dat <- bind_rows(seir_inc_dat1, seir_inc_dat2,seir_inc_dat_comb) %>% mutate(virus=factor(virus,levels=variant_levels))
    
    
    p_inc <- ggplot(seir_inc_dat) + 
      geom_line(aes(x=t,y=inc,col=virus))+
      geom_vline(xintercept=c(180),linetype="dashed",col="#D55E00") +
      geom_vline(xintercept=c(0),linetype="dashed",col="#0072B2") +
      variant_color_scale +
      xlab("Time") +
      scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
      #scale_y_continuous(expand=c(0,0)) +
      ylab("Per capita incidence") +
      theme_overall + 
      theme(legend.position=c(0.8,0.8)) +
      theme_nice_axes+ theme_no_x_axis
  }
  
  
  return(list(p_inc=p_inc,p_compartments=p_compartments,seir_res=seir_sim,
              virus1_inc=virus1_inc,virus2_inc=virus2_inc,
              virus1_inc_repeats=virus1_inc_repeats, virus2_inc_repeats=virus2_inc_repeats,
              virus1_inc_determ=virus1_inc_determ,
              virus2_inc_determ=virus2_inc_determ,
              virus1_inc_repeats_determ=virus1_inc_repeats_determ,
              virus2_inc_repeats_determ=virus2_inc_repeats_determ))
}


simulate_odin_2strain_seir <- function(sir_filename, tmax, pars, n_threads=4,seed=1,n_repeats=10){
  gen_sir <- odin.dust::odin_dust(sir_filename)
  
  dt <- pars["dt"]
  
  n_steps <- tmax/dt
  
  k_E <- pars["k_E"] 
  k_I <- pars["k_I"] 
  S_ini <- pars["S_ini"]
  beta1 <- pars["beta1"]
  beta2 <- pars["beta2"]
  beta3 <- pars["beta3"]
  beta4 <- pars["beta4"]
  tdur1 <- pars["tdur1"]
  tdur2 <- pars["tdur2"]
  tdur3 <- pars["tdur3"]
  tdur4 <- pars["tdur4"]
  transition <- pars["beta_transition_dur"]
  
  n_strains <- 4
  gamma <- pars["gamma"]
  sigma <- pars["sigma"]
  wane_rate <- pars["wane_rate"]
  
  seed_time1 <- pars["importtime1"]/dt
  seed_time2 <- pars["importtime2"]/dt
  
  seed_value1 <- rep(pars["seed_size1"],pars["seed_dur1"])
  seed_value2 <- rep(pars["seed_size2"],pars["seed_dur2"])
  
  strain1_trans <- pars["strain1_trans"]
  strain2_trans <- pars["strain2_trans"]
  strain12_trans <- pars["strain12_trans"]
  strain21_trans <- pars["strain21_trans"]
  
  strain12_immunity <- pars["strain12_immunity"]
  strain21_immunity <- pars["strain21_immunity"]
  
  beta_step <- c(rep(beta1,tdur1/dt),
                 seq(beta1, beta2, length.out=transition/dt),
                 rep(beta2,tdur2/dt),
                 seq(beta2, beta3, length.out=transition/dt),
                 rep(beta3,tdur3/dt),
                 seq(beta3, beta4, length.out=transition/dt),
                 rep(beta4,tdur4/dt))
  #beta_step <- rep(0.18, 10000)
  sir_model <- gen_sir$new(pars = list(dt = dt,
                                       S_ini = S_ini,
                                       beta_step = beta_step,
                                       n_strains = n_strains,
                                       gamma = gamma,
                                       sigma=sigma,
                                       waning_rate=wane_rate,
                                       trans_increase=rep(1, k_I),
                                       k_E=k_E,
                                       k_I=k_I,
                                       seed_step_start=seed_time1,
                                       seed_value=seed_value1,
                                       strain_seed_step_start=seed_time2,
                                       strain_seed_value=seed_value2,
                                       strain_transmission=c(strain1_trans,strain2_trans,strain12_trans,strain21_trans),
                                       cross_immunity=c(strain12_immunity,strain21_immunity)),
                           step = 1,
                           n_particles = n_repeats,
                           n_threads = n_threads,
                           seed = seed)
  
  n_particles <- n_repeats
  x <- array(NA, dim = c(sir_model$info()$len, n_particles, n_steps))
  for (t in seq_len(n_steps)) {
    x[ , , t] <- sir_model$run(t)
  }
  
  res <- plot_odin_output(x, sir_model)
  res$plot <- res$plot +
    theme_overall +
    theme_nice_axes+ theme_no_x_axis +
    theme(legend.position=c(0.8,0.8)) 
    #geom_line(data=tibble(x=seq(0,tdur1+tdur2+tdur3+tdur4 + (transition*3) - dt,by=dt),y=beta_step*100000),aes(x=x,y=y)) + scale_x_continuous(limits=c(0,tmax))
  res
}