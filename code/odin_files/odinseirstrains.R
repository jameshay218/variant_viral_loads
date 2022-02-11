## i for strain
## j for progression
##########################################
## SETUP
##########################################
## Default from odin example
## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Get total pop size
N <- S + sum(E) + sum(I) + sum(R)

##########################################
## SEEDING
##########################################
## Seeding of first variant: this will happen on the S->E flow
seed_step_start <- user()
seed_value[] <- user()
dim(seed_value) <- user()

seed_step_end <- seed_step_start + length(seed_value)
seed_rate <- if (step >= seed_step_start && step < seed_step_end)
    seed_value[as.integer(step - seed_step_start + 1)] else 0
seed <- rpois(seed_rate)

## New variant/strain seeding
strain_seed_step_start <- user()
strain_seed_value[] <- user()
dim(strain_seed_value) <- user()

strain_seed_step_end <- strain_seed_step_start + length(strain_seed_value)
strain_seed_rate <-
    if (step >= strain_seed_step_start && step < strain_seed_step_end)
        strain_seed_value[as.integer(step - strain_seed_step_start + 1)] else 0
strain_seed <- rpois(strain_seed_rate)



##########################################
## CORE EQUATIONS
##########################################
## Core equations for transitions between compartments:
update(S) <- S - out_S + in_S ## People leave S to become E of some variant
update(E[,]) <- E[i,j] + delta_E[i,j]
update(I[,]) <- I[i,j] + delta_I[i,j]
update(R[]) <- R[i] + delta_R[i]

## Cumulative incidence is those leaving S into E
update(obs_inc[]) <- inc[i]
update(obs_lambda[]) <- (S*(1 - exp(-(lambda[i]) * dt))/S_ini) + R[i]*p_R_progress[i]*(1-p_RS[i])/S_ini

##########################################
## TRANSITION PROBS
##########################################
## Individual probabilities of transition:
p_SE <- 1 - exp(-sum(lambda[]) * dt) # S to E based on FOI from all strains
p_E_progress <- 1 - exp(-sigma * k_E * dt) # progression of latent period
p_I_progress <- 1 - exp(-gamma * k_I * dt) # progression of infectious period

## Dimensions to account for strains
#dim(p_E_progress) <- c(n_strains, k_E)
#dim(p_I_progress) <- c(n_strains, k_I)
dim(n_S_progress) <- n_strains

## Draws from binomial distributions for numbers changing between compartments:
## Random number of individuals become exposed based on FOI lambda
## No one can move from S to E3 or E4 ie. can't go straight from S to being superinfected
n_S_progress_tot <- rbinom(S, p_SE) 
n_S_progress[] <- if(i == 1 || n_real_strains == 1) rbinom(n_S_progress_tot, rel_foi_strain[i]) else (if(i==2) n_S_progress_tot - n_S_progress[1] else 0) 

n_S_progress[1] <- n_S_progress[1] + min(S - n_S_progress_tot, seed) 
## Add seeds and make sure we don't remove all susceptibles
n_S_progress[2:n_strains] <-
    if (i < 3) min(n_S_progress[i] + strain_seed,
                   n_S_progress[i] + S -
                       n_S_progress_tot) else 0

out_S <- sum(n_S_progress)

##########################################
## SIMULATED TRANSITION NUMBERS
##########################################
n_E_progress[,] <- rbinom(E[i,j], p_E_progress) ## This is the number of people moving through E
n_I_progress[,] <- rbinom(I[i,j], p_I_progress) ## Number moving through I

## RECOVEREDS
## People go S->E->I->R, but can have "superinfection" by going R->E with respect
## to a particular strain. So can go Recovered with strain 1 to Exposed with strain 2.
## So people leave R compartment either by having waning immunity (will probably leave
## this as 0), or by going into one of the E superinfection compartments
## rate of progressing from R w/o superinfection is just waning_rate
## E1 is exposed strain 1 only
## E2 is exposed strain 2 only
## E3 is recovered strain 1, exposed strain 2
## E4 is recovered strain 2, exposed strain 1
## Note that (if n_real_strains == 2)
## cross_immunity[1] is the cross immunity of strain 1 against strain 2
## cross_immunity[2] is the cross immunity of strain 2 against strain 1
rate_R_progress[] <- waning_rate +
   if (n_strains == 1 || i > 2) 0 else
        lambda[3 - i] * (1 - cross_immunity[i])

p_R_progress[] <- 1 - exp(-rate_R_progress[i] * dt)
n_R_progress[] <- rbinom(R[i], p_R_progress[i])

## Waning
## Number becoming susceptible is proportion leaving R3/4 from waning
p_RS[] <- if (n_strains == 1 || i > 2) 1 else
    waning_rate/(waning_rate + lambda[3 - i] * (1 - cross_immunity[i]))
n_RS[] <- rbinom(n_R_progress[i], p_RS[i])

## New number in S is just all those who were waning
in_S <- sum(n_RS)

## Superinfection -- this is just those leaving R who weren't due to waning
n_RE[] <- n_R_progress[i] - n_RS[i]

delta_R[] <- n_I_progress[i,k_I]  - n_R_progress[i]

##########################################
## SORT OUT PRIMARY AND REINFECTIONS
##########################################
## For first E class, just number leaving S plus number of superinfections
## Otherwise, just stepping through
aux_E[,] <- (if (j == 1) n_S_progress[i] +
                 (if (i > 2) n_RE[i-2] else 0)
             else n_E_progress[i,j - 1]) - n_E_progress[i, j]
delta_E[,] <- aux_E[i,j]

inc[] <- n_S_progress[i] +
    (if (i > 2) n_RE[i-2] else 0)
dim(inc) <- n_strains

## INFECTED
## This is the number of people who move from each I stage
aux_I[,] <- (if (j == 1) n_E_progress[i, k_E]  ## Number moving out of the last E class
            else n_I_progress[i,j - 1]) - n_I_progress[i,j] ## Or the number moving into the class, less the number moving out
delta_I[,] <- aux_I[i,j]


## Now calculate the force of infection, which is the sum of contributions of all infected individuals modified by their strain-based modifier
lambda[] <- beta * strain_transmission[i] * sum(I[i,])/N


## Relative force of infection
rel_foi_strain[] <- (if(sum(lambda[]) == 0)
    (if(i==1) 1 else 0) else 
        min(lambda[i]/sum(lambda[]), as.numeric(1)))
dim(rel_foi_strain) <- n_real_strains


## Initial states are all zeroed as we will provide a state vector
## setting S and I based on the seeding model.
initial(S) <- S_ini
initial(E[,]) <- 0
#initial(E[1,1]) <- 100
initial(I[,]) <- 0
initial(R[]) <- 0
initial(obs_inc[]) <- 0
initial(obs_lambda[]) <- 0

## Parameters
## User defined parameters - default in parentheses:
beta_step[] <- user()
dim(beta_step) <- user()
## What we really want is min(step + 1, length(beta_step)) but that's not
## supported by odin (it could be made to support this).
beta <- if (as.integer(step) >= length(beta_step))
    beta_step[length(beta_step)] else beta_step[step + 1]

## Useful for debugging
initial(beta_out) <- beta_step[1]
update(beta_out) <- beta

gamma <- user(0.1)
sigma <- user(0.1)
waning_rate <- user(0.0)

k_E <- user()
k_I <- user()
S_ini <- user()


## Strain transmission
## I think this means we have "4" strains to cover coinfections, but
## this means we actually only have 2 real strains.
## So strain_transmission allows us to modify reinfection transmissibility
strain_transmission[] <- user()
dim(strain_transmission) <- n_strains
n_strains <- user()
n_real_strains <- if(n_strains == 4) 2 else 1

## Dimensions
dim(E) <- c(n_strains, k_E)
dim(I) <- c(n_strains, k_I)
dim(obs_inc) <- n_strains
dim(obs_lambda) <- n_real_strains


dim(n_E_progress) <- c(n_strains, k_E)
dim(aux_E) <- c(n_strains, k_E)
dim(delta_E) <- c(n_strains, k_E)

dim(n_I_progress) <- c(n_strains, k_E)
dim(aux_I) <- c(n_strains, k_E)
dim(delta_I) <- c(n_strains, k_E)

dim(R) <- n_strains
dim(delta_R) <- n_strains
dim(n_RS) <- n_strains
dim(p_RS) <- n_strains
dim(n_RE) <- n_strains
dim(p_R_progress) <- n_strains
dim(rate_R_progress) <- n_strains
dim(n_R_progress) <- n_strains

dim(cross_immunity) <- n_real_strains
cross_immunity[] <- user()

dim(lambda) <- n_real_strains