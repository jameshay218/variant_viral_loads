setwd("~/Documents/GitHub/driftSim/")
library(tidyverse)
library(getTBinR)
Rcpp::compileAttributes()
devtools::load_all()

set.seed(1)

## Read in individual-level parameters from previous model
pars <- read_csv("~/Google Drive/nCoV/sims_for_brian/swab_switch_SEIR_new/used_pars_swab_switch_SEIR_new_1.csv")

## Make peak viral load later
pars <- pars %>% dplyr::filter(!is.na(viral_peak))
pars <- pars[,c("tshift","mode","incu","t_wane","viral_peak")]
## t_wane<0.25 seems to mess things up
pars[,"t_wane"] <- pmax(0.25, pars[,"t_wane"])
## Make incubation period longer
pars <- as.matrix(pars)
## Symptomatic proportion is 30%
symptomatic_proportion <- 0.3
pars <- cbind(pars,rbernoulli(nrow(pars),symptomatic_proportion))

## Simulation control
flags <- c(1,1,1,1)
popN <- 1000000
## No deaths
hostPopn <- c(popN,0,5, 1/(10e10*365),1/7,0,10000)
n_variants <- 1
N <- popN*n_variants

## Infectiousness parameters
infectiousnessPars <- cbind(c(0.15),rep(1,n_variants),rep(7,n_variants))
crossImmunity <- matrix(c(1),nrow=n_variants)
seeds <- matrix(c(0,0,100),nrow=1,byrow=TRUE)

## Run simulation
run_simulation(flags, hostPopn,seeds,pars,infectiousnessPars,crossImmunity, start=0,end=100,VERBOSE=TRUE)

hosts <- read.csv("hosts.csv")
viruses <- read.csv("voutput.csv")

## Function to convert the saved host info into enumerate infections
parse_inf_hists <- function(hosts){
    all_hosts <- vector(mode = "list", length = nrow(hosts))
    index <- 1
    for(i in 1:nrow(hosts)){
        if(i %% 1000 == 0) print(i)
        inf_hist <- hosts$inf_history[i]
        inf_hist <- unlist(str_split(inf_hist, "-"))
        
        if(hosts$state[i] == 1) {
            inf_hist <- c(hosts$cur_inf[i], inf_hist)
        }
        
        for(j in 1:length(inf_hist)){
            inf_id <- as.numeric(inf_hist[j])
            if(!is.na(inf_id)){
                all_hosts[[index]] <- cbind(hosts[i,],"vid"=inf_id)
                index <- index + 1
            }
        }
        
    }
    all_hosts <- do.call("bind_rows",all_hosts)
    all_hosts
}

## Merge with virus objects to get full info
hosts_parsed <- parse_inf_hists(hosts)
comb <- hosts_parsed %>% full_join(viruses,by="vid") %>% as_tibble()

set.seed(1)
## Calculate offspring distribution over time and by age
with_offspring <- unique(comb$parentid)
no_offspring <- comb %>% dplyr::filter(!(vid %in% with_offspring)) %>% dplyr::select(vid,birth,age) %>% mutate(n=0) %>% rename(parentid=vid)
offspring <- comb %>% group_by(parentid) %>% tally() %>% left_join(comb %>% dplyr::select(vid,birth,age) %>% distinct() %>% rename(parentid=vid))
offspring <- bind_rows(offspring, no_offspring)
offspring <- offspring %>% dplyr::filter(birth > 0) %>% mutate(birth_week = floor(birth/7)) 
offspring_summaries <- offspring %>% group_by(birth_week,age) %>% summarize(y=mean(n),lower=quantile(n, 0.025),upper=quantile(n,0.975),N=n(),ysd=sd(n))

## SIMULATE REPORTING PROCESS
prob_report_early <- 0.1
prob_report_late <- 0.5
prob_report_switch <- 75
prob_report_weekend <- 0.6
t_report1 <- 77
t_report2 <- 100

## Reporting delay distribution
gamma_dist_mean_var_to_shape <- function(mean, var){
    scale <- var/mean
    shape <- mean/scale
    return(shape)
}
gamma_dist_mean_var_to_scale <- function(mean, var){
    scale <- var/mean
    shape <- mean/scale
    return(scale)
}
report_pars <- tibble(onset_time=0:150,mean1=7*(1 - 0:150*0.002),sd1=3*(1 - 0:150*0.002)) %>% 
    group_by(onset_time) %>%
    mutate(shape=gamma_dist_mean_var_to_shape(mean1,sd1),scale=gamma_dist_mean_var_to_scale(mean1,sd1))

to_process <- comb %>%
    mutate(onset_time=birth+ceiling(to)) %>% 
    ungroup() %>%
    left_join(report_pars) %>% 
    mutate(report_delay =extraDistr::rdgamma(n(), shape, scale=scale)) %>%
    mutate(report_time = onset_time + report_delay) %>%
    mutate(report_rate = ifelse(onset_time<=prob_report_switch,prob_report_early, prob_report_late)) %>%
    mutate(reported=rbernoulli(n(), report_rate)*(symptomatic.y==1)) %>%
    mutate(report_time=ifelse((report_time%% 7 %in% c(5,6))*rbernoulli(n(),1-prob_report_weekend),ceiling(report_time/7)*7,report_time)) %>%
    mutate(report_date = as.Date(report_time,origin="2022-07-01")) %>% 
    mutate(onset_date = as.Date(onset_time,origin="2022-07-01")) %>%
    left_join(age_pars %>% mutate(age=0:(n()-1))) %>%
    left_join(offspring %>% rename(vid=parentid)) %>% 
    mutate(infection_date = as.Date(birth,origin="2022-07-01"))

to_process$age_group <- factor(to_process$age_group,levels=age_pars$age_group)

######################################################
## PLOTS
######################################################
p_age_dist_1 <- ggplot(age_pars%>%mutate(i=0:(n()-1))) +
    geom_bar(data=to_process %>% dplyr::filter(reported==1) %>% dplyr::filter(report_date <= as.Date(t_report1-7,origin="2022-07-01"))  %>% group_by(age_group) %>% tally(),
              aes(x=age_group,y=n),stat="identity",col="black") + 
    #geom_line(aes(x=i+1,y=age_dist),col="blue",size=1.2) + 
    scale_y_continuous(expand=c(0,0)) +
    #scale_x_discrete(breaks=seq(0,8,by=1),labels=age_pars$age_group) +
    ylab("Count") +
    xlab("Age group") +
    scale_fill_who() +
    theme_who() +
    theme(axis.title.x=element_text(vjust=5,size=8),
          axis.text.x=element_text(size=7),axis.text.y=element_text(size=7),
          axis.title.y=element_text(size=8))
ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_age_dist1.png",p_age_dist_1, height=2.5,width=4.25,dpi=300,units="in")

p_age_dist_2 <- ggplot(age_pars%>%mutate(i=0:(n()-1))) +
    geom_bar(data=to_process %>%dplyr::filter(reported==1) %>% dplyr::filter(report_date <= as.Date(t_report2-7,origin="2022-07-01"))  %>% group_by(age_group) %>% tally(),
             aes(x=age_group,y=n),stat="identity",col="black") + 
    #geom_line(aes(x=i+1,y=age_dist),col="blue",size=1.2) + 
    scale_y_continuous(expand=c(0,0)) +
    #scale_x_discrete(breaks=seq(0,8,by=1),labels=age_pars$age_group) +
    ylab("Count") +
    xlab("Age group") +
    scale_fill_who() +
    theme_who() +
    theme(axis.title.x=element_text(vjust=5,size=8),
          axis.text.x=element_text(size=7),axis.text.y=element_text(size=7),
          axis.title.y=element_text(size=8))
ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_age_dist2.png",p_age_dist_2, height=2.5,width=4.25,dpi=300,units="in")

p_first_epi <- to_process %>%dplyr::filter(reported==1) %>% 
    dplyr::filter(report_date <= as.Date(t_report1,origin="2022-07-01")) %>%
    group_by(report_date, age, age_group) %>% 
    tally() %>% 
    rename(`Age group`=age_group) %>%
    ggplot() + 
    geom_bar(aes(x=report_date,y=n,fill=`Age group`),stat="identity") +
    #geom_line(data=viruses%>%group_by(birth)%>%tally(),aes(x=birth,y=n),col="red") + 
    scale_fill_who(palette="main") +
    theme_who() +
    scale_y_continuous(expand=c(0,0),limits=c(0,40),breaks=seq(0,40,by=10)) +
    #scale_x_date(limits=as.Date(c(0,77),origin="2022-07-01")) +
    ylab("Cases reported") +
    xlab("Date") +
    theme(legend.position = c(0.25,0.6),
          axis.title.x=element_text(vjust=5))+
    guides(fill=guide_legend(ncol=2))

epi_curve_1 <- to_process %>%dplyr::filter(reported==1) %>% 
    dplyr::filter(report_date <= as.Date(t_report1,origin="2022-07-01")) %>%
    group_by(report_date, age,age_group) %>% 
    rename(`Age group`=age_group) %>%
    tally() 

write_csv(epi_curve_1,"~/Google Drive/Teaching/Pandemic Exercise/outputs/epi_curve_1.csv")

ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_first_epi.png",p_first_epi, height=3,width=6,dpi=300,units="in")

p_second_epi <- to_process %>%dplyr::filter(reported==1) %>% 
    dplyr::filter(report_date <= as.Date(t_report2,origin="2022-07-01")) %>%
    group_by(report_date, age,age_group) %>% 
    rename(`Age group`=age_group) %>%
    tally() %>% 
    ggplot() + 
    geom_bar(aes(x=report_date,y=n,fill=`Age group`),stat="identity") +
    geom_vline(xintercept=as.Date(t_report1,origin="2022-07-01"),linetype="dotted",col="grey40",size=1) +
    scale_fill_who(palette="main") +
    theme_who() +
    scale_y_continuous(expand=c(0,0),limits=c(0,625),breaks=seq(0,600,by=50)) +
    #scale_x_date(breaks="1 weeks", limits=as.Date(c("2022-07-11","2022-09-19"))) +
    ylab("Cases reported") +
    xlab("Date") +
    theme(legend.position = c(0.25,0.6),
          axis.title.x=element_text(vjust=5))+
    guides(fill=guide_legend(ncol=2))


epi_curve_2 <- to_process %>%dplyr::filter(reported==1) %>% 
    dplyr::filter(report_date <= as.Date(t_report2,origin="2022-07-01")) %>%
    group_by(report_date, age,age_group) %>% 
    rename(`Age group`=age_group) %>%
    tally() 

write_csv(epi_curve_2,"~/Google Drive/Teaching/Pandemic Exercise/outputs/epi_curve_2.csv")

ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_second_epi.png",p_second_epi,height=3,width=6,dpi=300,units="in")


## Generate cluster of cases to investigate
## Assume we only "contact trace" back up to 5 days before symptom onset
## 1. Find an index case with symptom onset on some day
## 2. Find all offspring of this index case that have arisin up to 5 days before symptom onset or after.
## 3. Report only symptomatic cases or those with viral load > 5 (or something)
## First outbreak investigation
index_case <- to_process %>%
    dplyr::filter(onset_date >= "2022-09-01", onset_date <= "2022-09-07", reported==1) %>%  arrange(-n)

vl_lod1 <- 3 ## Limit of detection of test
contact_trace_past1 <- 14 ## Test all infections from contacts 7 days prior to detect of the index case

contract_tracing_study <- function(n_studies, vl_lod1, contact_trace_past1, index_case, to_process){
    linelists1 <- NULL
    for(index in 1:n_studies){
        index_case_tmp <- index_case[index,] %>% mutate(generation=1)
        index_vid <- index_case_tmp %>% pull(vid)
        report_time_tmp <- index_case[index,] %>% pull(report_time)
        report_date_tmp <- index_case[index,] %>% pull(report_date)
        
        ## Second generation of infections
        ## Find viruses caused by this index case, calculate viral load at the detection time and see if detected
        second_gen <- to_process %>%dplyr::filter(parentid == index_vid) %>% 
            mutate(report_time = report_time_tmp)%>% 
            mutate(report_date = report_date_tmp) %>% 
            mutate(t_since_infection=as.numeric(report_date - infection_date))%>%
            rowwise() %>% 
            mutate(vl=solve_viral_kinetics(c(0,t_since_infection), tg,tp,to,tw,alpha)[2]) %>%
            mutate(vl = pmin(vl, 20)) %>%
           dplyr::filter(symptomatic.y==1 | (symptomatic.y==0 & vl >= vl_lod1 & parent_age_at_creation <= contact_trace_past1)) %>% 
            mutate(generation = 2)
        second_gen_index <- second_gen %>% pull(vid)
        
        ## Third generation of infections
        third_gen <- to_process %>%dplyr::filter(parentid %in% second_gen_index) %>%
            mutate(report_time = report_time_tmp)%>% 
            mutate(report_date = report_date_tmp) %>% 
            mutate(t_since_infection=as.numeric(report_date - infection_date)) 
        
        ## If no third generation of infections, then pass
        if(nrow(third_gen) > 0){
            third_gen <- third_gen %>%
                rowwise() %>%
            mutate(vl=solve_viral_kinetics(c(0,t_since_infection), tg,tp,to,tw,alpha)[2]) %>%
            mutate(vl = pmin(vl, 20)) %>%
            mutate(generation=3) %>%
           dplyr::filter(symptomatic.y==1 | (symptomatic.y==0 & vl >= vl_lod1 & parent_age_at_creation <= contact_trace_past1)) 
            third_gen_index <- third_gen %>% pull(vid)
        } else {
            third_gen <- NULL
        }
        
        ## If there was a third generation of infections, look for the fourth generation
        if(!is.null(third_gen)){
            fourth_gen <- to_process %>%dplyr::filter(parentid %in% third_gen_index) %>% 
                mutate(report_time = report_time_tmp)%>% 
                mutate(report_date = report_date_tmp) %>% 
                mutate(t_since_infection=as.numeric(report_date - infection_date)) 
        } else {
            fourth_gen <- NULL
        }
        if(!is.null(fourth_gen) && nrow(fourth_gen) > 0){
            fourth_gen <- fourth_gen %>%
            rowwise() %>%
            mutate(vl=solve_viral_kinetics(c(0,t_since_infection), tg,tp,to,tw,alpha)[2]) %>%
            mutate(vl = pmin(vl, 20)) %>%
           dplyr::filter(symptomatic.y==1 | (symptomatic.y==0 & vl >= vl_lod1 & parent_age_at_creation <= contact_trace_past1))%>% 
            mutate(generation = 4)
        } else {
            fourth_gen <- NULL
        }
        
        ## Combine all detected infections
        linelists1[[index]] <- bind_rows(index_case_tmp, second_gen, third_gen,fourth_gen)  %>%
            mutate(outbreak=index) %>% 
            mutate(report_time = report_time_tmp)%>% 
            mutate(report_date = report_date_tmp)
    }
    linelist1 <- do.call("bind_rows",linelists1) %>% rename(R=n) %>% mutate(t_since_onset = report_time - onset_time) %>%
        dplyr::select(!c(vaccine_time, inf_history, variant_id, variant,infectionNo,tg,tp,tw,alpha,
                         infectiousnessMax,infectiousnessGradient,infectiousnessInflection, host_vacc_status,
                         age_dist,infectivity, susceptibility,state,last_vid,cur_inf,symptomatic.x,age,mean1,sd1,shape,scale))
}
linelist1 <- contract_tracing_study(50, vl_lod1, contact_trace_past1, index_case, to_process)
## Get epi pars from linelist
## Calculate offspring distribution over time and by age
with_offspring_ll <- unique(linelist1$parentid)
no_offspring_ll <- linelist1 %>%dplyr::filter(generation < 4) %>% dplyr::filter(!(vid %in% with_offspring)) %>% dplyr::select(vid,birth,age_group,generation) %>% mutate(n=0) %>% rename(parentid=vid)
offspring_ll <- linelist1 %>%dplyr::filter(generation > 1) %>% group_by(parentid) %>% tally() %>% left_join(linelist1 %>% dplyr::select(vid,birth,age_group,generation) %>% distinct() %>% rename(parentid=vid))
offspring_ll <- bind_rows(offspring_ll, no_offspring_ll)
#offspring_ll <- offspring_ll %>% dplyr::filter(birth > 0) %>% mutate(birth_week = floor(birth/7)) 
offspring_summaries_ll_1 <- offspring_ll %>% group_by(age_group) %>% summarize(R=mean(n),lower=quantile(n, 0.025),upper=quantile(n,0.975),N=n(),ysd=sd(n))

p_r_number1 <- offspring_ll %>% 
    ggplot() + 
    geom_histogram(aes(x=n),binwidth=1,col="black",fill="grey40") + 
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits=c(min(offspring_ll$n-0.5),max(offspring_ll$n+0.5)),breaks=seq(min(offspring_ll$n), max(offspring_ll$n),by=2)) +
    facet_wrap(~age_group,scales="free") + 
    ylab("Count") +
    xlab("Number of secondary cases caused by index case") +
    theme_who() +
    theme(strip.text=element_text(size=8))+
    theme(axis.title.x=element_text(vjust=5,size=8),
          axis.text.x=element_text(size=7),axis.text.y=element_text(size=7),
          axis.title.y=element_text(size=8)) + theme(panel.spacing = unit(0.1, "lines")) 

p_vl_plot1 <- linelist1 %>%dplyr::filter(vl > vl_lod1) %>% mutate(vl = round(vl, 1)) %>% 
    rename(Symptomatic=symptomatic.y) %>%
    mutate(Symptomatic=ifelse(Symptomatic==1,"Yes",'No')) %>%
    ggplot() + 
    geom_hline(yintercept=vl_lod1,linetype="dashed") +
    geom_vline(xintercept=mean(linelist1$to),linetype="dotted",col="red") +
    geom_point(aes(x=t_since_infection,y=vl,col=`Symptomatic`),size=0.5) +
    scale_y_continuous(limits=c(1.8,11),expand=c(0,0),breaks=seq(0,11,by=1)) +
    scale_x_continuous(limits=c(0, 20)) +
    scale_color_manual(values=c("Yes"="red","No"="blue")) +
    xlab("Days since exposure") +
    ylab("Viral load (log10 RNA copies/ml)") +
    theme_who()+
    theme(axis.title.x=element_text(vjust=5,size=8),
          axis.text.x=element_text(size=7),axis.text.y=element_text(size=7),
          axis.title.y=element_text(size=8), legend.position=c(0.8,0.8))

#linelist <- linelist %>% select(state, age, vid, birth, death, parentid, parent_age_at_creation, symptomatic.y, onset_time, report_time, reported, age_group, R, outbreak)


## First outbreak investigation
index_case <- to_process %>%
    dplyr::filter(report_date >= "2022-10-01", report_date <= "2022-10-07", reported==1) %>%  arrange(-n)

vl_lod2 <- 2
contact_trace_past2 <- 14
linelist2 <- contract_tracing_study(50, vl_lod2, contact_trace_past2, index_case, to_process)

## Get epi pars from linelist
## Calculate offspring distribution over time and by age
with_offspring_ll <- unique(linelist2$parentid)
no_offspring_ll <- linelist2 %>%dplyr::filter(generation < 4) %>% dplyr::filter(!(vid %in% with_offspring)) %>% dplyr::select(vid,birth,age_group,generation) %>% mutate(n=0) %>% rename(parentid=vid)
offspring_ll <- linelist2 %>%dplyr::filter(generation > 1) %>% group_by(parentid) %>% tally() %>% left_join(linelist2 %>% dplyr::select(vid,birth,age_group,generation) %>% distinct() %>% rename(parentid=vid))
offspring_ll <- bind_rows(offspring_ll, no_offspring_ll)
offspring_summaries_ll_2 <- offspring_ll %>% group_by(age_group) %>% summarize(R=mean(n),lower=quantile(n, 0.025),upper=quantile(n,0.975),N=n(),ysd=sd(n))

p_r_number2 <- offspring_ll %>% 
    ggplot() + 
    geom_histogram(aes(x=n),binwidth=1,col="black",fill="grey40") + 
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(limits=c(min(offspring_ll$n-0.5),max(offspring_ll$n+0.5)),breaks=seq(min(offspring_ll$n), max(offspring_ll$n),by=2)) +
    facet_wrap(~age_group,scales="free") + 
    ylab("Count") +
    xlab("Number of secondary cases caused by index case") +
    theme_who() +
    theme(strip.text=element_text(size=8))+
    theme(axis.title.x=element_text(vjust=5,size=8),
          axis.text.x=element_text(size=7),axis.text.y=element_text(size=7),
          axis.title.y=element_text(size=8)) + theme(panel.spacing = unit(0.1, "lines")) 

p_vl_plot2 <- linelist2 %>%dplyr::filter(vl > vl_lod2) %>% mutate(vl = round(vl, 1)) %>% 
    rename(Symptomatic=symptomatic.y) %>%
    mutate(Symptomatic=ifelse(Symptomatic==1,"Yes",'No')) %>%
    ggplot() + 
    geom_hline(yintercept=vl_lod2,linetype="dashed") +
    geom_vline(xintercept=mean(linelist1$to),linetype="dotted",col="red") +
    geom_point(aes(x=t_since_infection,y=vl,col=Symptomatic),size=0.5) +
    scale_y_continuous(limits=c(1.8,11),expand=c(0,0),breaks=seq(0,11,by=1)) +
    scale_x_continuous(limits=c(0, 20)) +
    scale_color_manual(values=c("Yes"="red","No"="blue")) +
    xlab("Days since exposure") +
    ylab("Viral load (log10 RNA copies/ml)") +
    theme_who()+
    theme(axis.title.x=element_text(vjust=5,size=8),
          axis.text.x=element_text(size=7),axis.text.y=element_text(size=7),
          axis.title.y=element_text(size=8), legend.position=c(0.8,0.8))


ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_r_dist1.png",p_r_number1,height=3.5,width=5,dpi=300,units="in")
ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_r_dist2.png",p_r_number2,height=3.5,width=5,dpi=300,units="in")
ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_vl_plot1.png",p_vl_plot1,height=2.5,width=4.5,dpi=300,units="in")
ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_vl_plot2.png",p_vl_plot2,height=2.5,width=4.5,dpi=300,units="in")


generate_symptoms <- function(mean_n_symptoms, symptom_list,symptom_probs){
    symptoms <- sample(symptom_list, max(1,min(rpois(1,mean_n_symptoms), length(symptom_list))),replace=FALSE,prob=symptom_probs)
    if(length(symptoms) > 0){
        symptoms <- paste0(symptoms, collapse=", ")
    } else {
        symptoms <- "none"
    }
    return(symptoms)
}
symptom_list <- c("nausea","bloating","diarrhea","extreme fatique","fever","vomiting","headache","body aches")
symptom_probs <- c(0.3,0.4,0.4,0.8,0.1,0.25,0.05,0.1)

linelist_report1 <- linelist1 %>% select(infection_date, symptomatic.y,to,onset_date, report_date,age_group,R, vl,t_since_infection) %>% 
    rowwise() %>% 
    mutate(symptoms=ifelse(symptomatic.y==1, generate_symptoms(3, symptom_list, symptom_probs), "None")) %>%
    mutate(to = ifelse(symptomatic.y == 1 & !is.na(t_since_infection),round(to,0), NA)) %>%
    select(!c(symptomatic.y)) %>%
    mutate(vl=round(vl, 1)) %>%
    arrange(infection_date) %>%
    rename(`Infection date`=infection_date, `Pre-symptomatic duration (days)`=to, `Symptom onset date`=onset_date,
           `Report date`=report_date, `Age group`=age_group,`Secondary cases identified`=R,`Viral load`=vl,`Time since infection (days)`=t_since_infection,
           `Symptoms`=symptoms)

write_csv(linelist_report1, "~/Google Drive/Teaching/Pandemic Exercise/outputs/linelist_report_1.csv")

linelist_report2 <- linelist2 %>% select(infection_date, symptomatic.y,to,onset_date, report_date,age_group,R, vl,t_since_infection) %>% 
    rowwise() %>% 
    mutate(symptoms=ifelse(symptomatic.y==1, generate_symptoms(3, symptom_list, symptom_probs), "none")) %>%
    mutate(to = ifelse(symptomatic.y == 1 & !is.na(t_since_infection),round(to,0), NA)) %>%
    select(!c(symptomatic.y)) %>%
    mutate(vl=round(vl, 1)) %>%
    arrange(infection_date) %>%
    rename(`Infection date`=infection_date, `Pre-symptomatic duration (days)`=to, `Symptom onset date`=onset_date,
           `Report date`=report_date, `Age group`=age_group,`Secondary cases identified`=R,`Viral load`=vl,`Time since infection (days)`=t_since_infection,
           `Symptoms`=symptoms)

write_csv(linelist_report2, "~/Google Drive/Teaching/Pandemic Exercise/outputs/linelist_report_2.csv")

############################
## SIMULATION TRUTH PLOTS
############################
p_full_epi_curve <- to_process %>%dplyr::filter(reported==1) %>% 
    dplyr::filter(report_date <= as.Date(100,origin="2022-07-01")) %>%
    group_by(report_date, age,age_group) %>% 
    rename(`Age group`=age_group) %>%
    tally() %>% 
    ggplot() + 
    geom_bar(aes(x=report_date,y=n,fill=`Age group`),stat="identity") +
    scale_fill_who(palette="main") +
    scale_y_continuous(expand=c(0,0),limits=c(0,625),breaks=seq(0,600,by=50)) +
    xlab("Date") +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x=element_text(vjust=5))+
    guides(fill=guide_legend(ncol=2)) +
    geom_line(data=to_process %>% group_by(infection_date) %>% tally(), aes(x=infection_date,y=n),col="red",size=0.75) + 
    geom_line(data=to_process %>%dplyr::filter(symptomatic.y==1) %>% group_by(infection_date) %>% tally(), aes(x=infection_date,y=n),col="blue",size=0.75) + 
    scale_y_continuous(limits=c(0,13000),expand=c(0,0),breaks=seq(0,13000,by=2000)) +
    geom_hline(yintercept=600,size=0.25,col="grey40") +
    ylab("True incidence of infections")

epi_curve_full <- to_process %>%
    dplyr::filter(onset_date <= as.Date(t_report2,origin="2022-07-01")) %>%
    group_by(onset_date) %>% 
    tally() 

epi_curve_full_infections <- to_process %>%
    dplyr::filter(infection_date <= as.Date(t_report2,origin="2022-07-01")) %>%
    group_by(infection_date) %>% 
    tally() 

write_csv(epi_curve_full,"~/Google Drive/Teaching/Pandemic Exercise/outputs/full_epi_curve.csv")
write_csv(epi_curve_full_infections,"~/Google Drive/Teaching/Pandemic Exercise/outputs/full_epi_curve_infections.csv")

ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_full_epi.png",p_full_epi_curve,height=3,width=6,dpi=300,units="in")


age_pars$age_group <- factor(age_pars$age_group,levels=age_pars$age_group)
p_age_infectivity <- ggplot(age_pars) +
    geom_bar(aes(x=age_group,y=infectivity),stat="identity") +
    theme_bw() +
    xlab("Age group") +
    ylab("Infectiousness relative to maximum") +
    scale_y_continuous(expand=c(0,0))
p_age_susceptibility <- ggplot(age_pars) +
    geom_bar(aes(x=age_group,y=susceptibility),stat="identity") +
    theme_bw() +
    xlab("Age group") +
    ylab("Susceptibility relative to maximum") +
    scale_y_continuous(expand=c(0,0))


ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_age_infect.png",p_age_infectivity,height=4,width=6,dpi=300,units="in")
ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_age_suscept.png",p_age_susceptibility,height=4,width=6,dpi=300,units="in")


true_vl <- t(apply(pars[1:1000,], 1, function(x) solve_viral_kinetics(seq(0,35,by=1),x[1],x[2],x[3],x[4],x[5])))
true_vl <- true_vl + rnorm(length(true_vl),0,pars[,6])
true_vl[true_vl < 0] <- 0

true_vl <- reshape2::melt(true_vl)
colnames(true_vl) <- c("i","t","vl")
p_true_vl <- ggplot(true_vl) + geom_point(aes(x=t,y=vl),size=0.25) + scale_x_continuous(breaks=seq(0,35,by=5))+ scale_y_continuous(expand=c(0,0),breaks=seq(0,12,by=2),limits=c(0,12)) +
    ylab("Viral load") + xlab("Days post infection") + theme_bw()

dat_prop_detect <- true_vl %>% mutate(detectable=vl>=3) %>% drop_na() %>% group_by(t)%>% summarize(n_detect=sum(detectable),N=n(),prop=n_detect/N)
p_true_detect <- ggplot(dat_prop_detect) + geom_line(aes(x=t,y=prop)) + scale_x_continuous(breaks=seq(0,35,by=5))+ 
    ylab("Proportion of infections detectable") + xlab("Days post infection") + theme_bw()


infectivity_profiles <- t(apply(pars[1:25,], 1, function(x) 5*solve_infectiousness(seq(0,20,by=1),x[1],x[2],x[3],x[4],x[5],0.28,1,7)))
infectivity_profiles <- reshape2::melt(infectivity_profiles)
colnames(infectivity_profiles) <- c("i","t","infectiousness")
p_true_infectivity <- ggplot(infectivity_profiles) + geom_line(aes(x=t,y=infectiousness,group=i),size=0.25,col="grey50") + 
    scale_x_continuous(breaks=seq(0,35,by=5))+ 
    scale_y_continuous(expand=c(0,0),limits=c(0,1.5)) +
    ylab("Expected number of infected contacts") + xlab("Days post infection") + theme_bw() + ggtitle("Infectiousness in oldest age group")

p_true_infectivity_young <- ggplot(infectivity_profiles) + geom_line(aes(x=t,y=infectiousness*0.3,group=i),size=0.25,col="darkgreen") + scale_x_continuous(breaks=seq(0,35,by=5))+ 
    scale_y_continuous(expand=c(0,0),limits=c(0,1.5)) +
    ylab("Expected number of infected contacts") + xlab("Days post infection") + theme_bw() + ggtitle("Infectiousness in 20-29 age group")

p_r_numbers_true <- ggplot(offspring_summaries %>% 
           left_join(age_pars %>% mutate(age=0:(n()-1))) %>%
           left_join(offspring %>% mutate(birth=as.Date(birth,origin="2022-07-01")) %>% group_by(birth_week) %>% summarize(mean_week=mean(birth)))) + 
    geom_hline(yintercept=1,linetype="dashed") +
    geom_line(aes(x=mean_week,y=y,col=age_group)) + 
    geom_vline(xintercept=as.Date("2022-09-25")) +
    ylab("Effective reproductive number") +
    xlab("Date") +
    facet_wrap(~age_group) +
    theme_bw() +
    scale_color_viridis_d() +
    theme(legend.position = "none")

p_incu_periods <- ggplot(to_process) +
    geom_histogram(aes(x=to,y=..density..),binwidth=1,col="black") +
    theme_bw() +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(breaks=seq(0,50,by=5)) +
    xlab("Delay between infection and symptom onset") +
    ylab("Proportion")

p_reporting_proportion <- to_process %>% select(onset_date,report_rate) %>% distinct() %>%
    ggplot() + geom_line(aes(x=onset_date,y=report_rate)) +
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1),expand=c(0,0)) + theme_bw() +
    xlab("Symptom onset date") +
    ylab("Probability of being reported")

p_reporting_delay <- report_pars %>% mutate(onset_date=as.Date(onset_time,origin="2022-07-01")) %>%
    ggplot() + 
    geom_ribbon(aes(x=onset_date,ymin=mean1-sd1,ymax=mean1+sd1),fill="blue",alpha=0.25) +
    geom_line(aes(x=onset_date,y=mean1)) +
    theme_bw() +
    scale_y_continuous(limits=c(0,10),breaks=seq(0,10,by=2)) +
    xlab("Symptom onset date") +
    ylab("Delay between symptom onset\n and being tested")


ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_true_detect.png",p_true_detect,height=3,width=6,dpi=300,units="in")
ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_true_vl.png",p_true_vl,height=3,width=6,dpi=300,units="in")
ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_true_infectivity.png",p_true_infectivity,height=3,width=6,dpi=300,units="in")
ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_true_infectivity_young.png",p_true_infectivity_young,height=3,width=6,dpi=300,units="in")
ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_true_r_numbers.png",p_r_numbers_true,height=5,width=6,dpi=300,units="in")
ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_true_incu_period.png",p_incu_periods,height=3,width=6,dpi=300,units="in")
ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_reporting_proportion.png",p_reporting_proportion,height=3,width=6,dpi=300,units="in")
ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_reporting_delay.png",p_reporting_delay,height=3,width=6,dpi=300,units="in")

## Generation interval distribution
p_generation_interval <- to_process %>%dplyr::filter(parent_age_at_creation > 0) %>% 
    ggplot() + 
    geom_histogram(aes(x=parent_age_at_creation),binwidth=1,col="black")+
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    ylab("Count") +
    xlab("Delay between infection in infector\n and infection in infectee")

ggsave("~/Google Drive/Teaching/Pandemic Exercise/figures/p_generation_interval.png",p_generation_interval,height=3,width=6,dpi=300,units="in")

## Save full linelist
to_save <- to_process %>% 
    mutate(recovery_date=as.Date(death,origin="2022-07-01")) %>%
    select(infection_date,death, vid, parentid,to,symptomatic.y, onset_date, report_delay, reported, report_date, age_group) %>%
    rename(`Infection date`=infection_date,`Recovery date`=death, `Virus ID`=vid, `Virus ID of infector`=parentid, 
           `Incubation period`=to, `Symptomatic?`=symptomatic.y,`Onset date`=onset_date,
           `Reporting delay`=report_delay,`Is reported?`=reported,`Report date`=report_date,
           `Age group`=age_group)
write_csv(to_save,"~/Google Drive/Teaching/Pandemic Exercise/outputs/full_linelist.csv")
