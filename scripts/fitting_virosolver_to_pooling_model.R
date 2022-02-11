
likelihood_test <- function(obs_dat, pars){
    LOD <- pars["intercept"]
    
    ## Because we have a different standard deviation for different times
    ## Time at which standard deviation is reduced
    t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
    sd_mod <- rep(pars["sd_mod"], length(ages)) ## Up until t_switch, full standard deviation
    
    ## Prior to t_switch, 1
    unmod_vec <- 1:min(t_switch,length(ages))
    sd_mod[unmod_vec] <- 1
    
    ## For the next sd_mod_wane days, decrease linearly
    decrease_vec <- (t_switch+1):(t_switch+pars["sd_mod_wane"])
    ## The rest are at sd_mod
    sd_mod[decrease_vec] <- 1 - ((1-pars["sd_mod"])/pars["sd_mod_wane"])*seq_len(pars["sd_mod_wane"])
    
    ## Get the modal Ct value
    viral_loads <- virosolver::viral_load_func(pars, ages)
    
    prob_detectable_dat <- sapply(ages[ages > 0], function(a) virosolver::prop_detectable_cpp(a, viral_loads[a],
                                                                                              pars["obs_sd"]*sd_mod[a], pars["intercept"],
                                                                                              t_switch, pars["prob_detect"]))
    
    detectable_prob <- function(X_i, a){
        virosolver::p_a(X_i, a, pars, viral_loads,sd_mod)*prob_detectable_dat[a]
    }
    undetectable_prob <- function(a) {
        1 - prob_detectable_dat[a]
    }
    
    lik_func <- function(ct, a){
        log((ct >= LOD)*undetectable_prob(a) + (ct < LOD)*detectable_prob(ct,a))
    }
    obs_dat %>% mutate(lik = lik_func(ct, t)) %>% pull(lik) -> lik
    sum(lik)
}
create_posterior_func <- function(parTab, data, ages=1:50,PRIOR_FUNC = NULL,func_version="posterior", ADD_NOISE=FALSE){
    par_names <- parTab$names
    pars <- parTab$values
    names(pars) <- par_names
    
    if(func_version == "posterior"){
        f <- function(pars){
            names(pars) <- par_names
            lik <- likelihood_test(data,pars)
            if(!is.null(PRIOR_FUNC)){
                lik <- lik + PRIOR_FUNC(pars)
            }
            lik
        }
    } else {
        f <- function(pars){
            names(pars) <- par_names
            LOD <- pars["intercept"]
            
            ## Because we have a different standard deviation for different times
            ## Time at which standard deviation is reduced
            t_switch <-  pars["t_switch"] + pars["desired_mode"] + pars["tshift"]
            sd_mod <- rep(pars["sd_mod"], length(ages)) ## Up until t_switch, full standard deviation
            
            ## Prior to t_switch, 1
            unmod_vec <- 1:min(t_switch,length(ages))
            sd_mod[unmod_vec] <- 1
            
            ## For the next sd_mod_wane days, decrease linearly
            decrease_vec <- (t_switch+1):(t_switch+pars["sd_mod_wane"])
            ## The rest are at sd_mod
            sd_mod[decrease_vec] <- 1 - ((1-pars["sd_mod"])/pars["sd_mod_wane"])*seq_len(pars["sd_mod_wane"])
            
            ## Get the modal Ct value
            viral_loads <- virosolver::viral_load_func(pars, ages)
            
            if(ADD_NOISE){
                viral_loads <- extraDistr::rgumbel(length(viral_loads),viral_loads, sd_mod*pars["obs_sd"])
            }
            viral_loads
        }
    }
    f
}

ages <- 1:50
x <- melted_cts %>% dplyr::filter(i <= 100, t > 0)
tmp <- likelihood_test(x, vl_pars)
microbenchmark::microbenchmark(likelihood_test(x, vl_pars),times=100)

parTab <- read.csv("~/Documents/GitHub/variant_viral_loads/pars/partab_vl_prior_fit.csv")
f <- create_posterior_func(parTab, x, NULL)
f(parTab$values)

library(lazymcmc)
mcmcPars <- c("iterations"=50000,"popt"=0.44,"opt_freq"=500,
                    "thin"=5,"adaptive_period"=2000,"save_block"=1000)

run_MCMC(parTab, data=x,mcmcPars=mcmcPars,CREATE_POSTERIOR_FUNC = create_posterior_func,filename="tmp")
chain <- read.csv("tmp_univariate_chain.csv")
chain <- chain[chain$sampno > mcmcPars["adaptive_period"],]
#plot(coda::as.mcmc(chain))


f <- create_posterior_func(parTab, x, PRIOR_FUNC=NULL,ages=1:50,func_version = "model", ADD_NOISE=TRUE)
plot(f(parTab$values))

y <- t(apply(chain[,2:ncol(chain)], 1, function(x) f(x)))
y[y > 40] <- 40
quants <- as_tibble(t(apply(y, 2, function(x) quantile(x[x < 40], c(0.025,0.5,0.975)))))
colnames(quants) <- c("lower","median","upper")
quants <- quants %>% mutate(t=1:50)
ggplot(quants) + geom_line(aes(x=t,y=median))

best_pars <- get_best_pars(chain)
f1 <- create_posterior_func(parTab, x, PRIOR_FUNC=NULL,ages=1:50,func_version = "model", ADD_NOISE=FALSE)
best_vl <- tibble(x=1:50,y=f1(best_pars))

ggplot() + 
    geom_line(data=quants,aes(x=t,y=median)) +
    geom_ribbon(data=quants,aes(x=t,ymin=lower,ymax=upper),alpha=0.25) +
    geom_line(data=best_vl,aes(x=x,y=y),col="blue") +
    geom_jitter(data=melted_cts %>% dplyr::filter(ct < 40, i <= 50), aes(x=t,y=ct),size=0.1,height=0.1,width=0.25) + 
    geom_line(data=melted_cts %>% dplyr::filter(ct < 40) %>% group_by(t) %>% summarize(mean_ct=median(ct)),aes(x=t,y=mean_ct),col="red") +
    scale_y_continuous(trans="reverse") + 
    #geom_line(data=omg,aes(x=x,y=y)) +
    scale_x_continuous(breaks=seq(0,50,by=5))

parTab1 <- parTab
parTab1$values <- best_pars

write.csv(parTab1,"~/Documents/GitHub/variant_viral_loads/pars/partab_vl_prior_FITTED.csv")
