## Takes a vector of model parameters, initial states and solve times and solves the two-strain SEIR model
run_2strain_seir_simulation <- function(pars, states, times){
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
  
  seir_inc_dat1 <- tibble(t=times,inc=virus1_inc,virus="Original variant")
  seir_inc_dat2 <- tibble(t=times,inc=virus2_inc,virus="New variant")
  seir_inc_dat_comb <- tibble(t=times,inc=virus1_inc+virus2_inc,virus="Overall")
  seir_inc_dat <- bind_rows(seir_inc_dat1, seir_inc_dat2,seir_inc_dat_comb) %>% mutate(virus=factor(virus,levels=variant_levels))
  
  
  p_inc <- ggplot(seir_inc_dat) + 
    geom_line(aes(x=t,y=inc,col=virus))+
    geom_vline(xintercept=c(180),linetype="dashed",col="red") +
    geom_vline(xintercept=c(0),linetype="dashed",col="blue") +
    variant_color_scale +
    xlab("Time") +
    scale_x_continuous(limits=c(0,max(times)),breaks=seq(0,550,by=50)) +
    #scale_y_continuous(expand=c(0,0)) +
    ylab("Per capita incidence") +
    theme_overall + 
    theme(legend.position=c(0.8,0.8)) +
    theme_nice_axes+ theme_no_x_axis
  return(list(p_inc=p_inc,p_compartments=p_compartments,seir_res=seir_sim,virus1_inc=virus1_inc,virus2_inc=virus2_inc))
}

