#gen_sir <- odin.dust::odin_dust("../code/odin_files/odinseirstrains.R")


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
    tdur1 <- pars["tdur1"]
    tdur2 <- pars["tdur2"]
    tdur3 <- pars["tdur3"]
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
                   rep(beta3,tdur3/dt))
    
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
    res
}

#res$dat %>% filter(time == 200, compartment=="obs_inc", strain=="1")


## Reporting delay distribution
# gamma_dist_mean_var_to_shape <- function(mean, var){
#     scale <- var/mean
#     shape <- mean/scale
#     return(shape)
# }
# gamma_dist_mean_var_to_scale <- function(mean, var){
#     scale <- var/mean
#     shape <- mean/scale
#     return(scale)
# }
#par(mfrow=c(2,1))
#plot(dgamma(0:25, shape=k_E,rate=k_E/3),type='l')
#plot(dgamma(0:25, shape=k_I,rate=k_I/7),type='l')
#par(mfrow=c(1,1))
#mean(rgamma(10000, shape=k_I,rate=k_I/7))
#pgamma(0,shape=k_I,scale=k_I*(1/7),lower.tail=FALSE)

