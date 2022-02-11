## Default from odin example
## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Seeding of first variant: this will happen on the S->E flow
#seed_step_end_var1 <- seed_start_var1 + seed_duration_var1
#seed_var1 <- if (step >= seed_start_var1 && step < seed_step_end_var1)
#    seed_value_var1 else 0

## Get user inputs for these
#seed_start_var1 <- user()
#seed_duration_var1 <- user()
seed_value_var1 <- user()
seed_var1 <- 0

N <- S + sum(E) + sum(I) + R

## Core equations for transitions between compartments:
update(S) <- S - delta_S ## People leave S to become E of some variant
update(E[]) <- E[i] + delta_E[i]
update(I[]) <- I[i] + delta_I[i]
update(R) <- R + delta_R
update(cumu_inc) <- delta_S

## Individual probabilities of transition:
p_S_progress <- 1 - exp(-lambda * dt) # S to I - age dependent
p_E_progress <- 1 - exp(-sigma * k_E * dt) # progression of latent period
p_I_progress <- 1 - exp(-gamma * k_I * dt) # progression of infectious period

## Draws from binomial distributions for numbers changing between
## compartments:
n_S_progress <- rbinom(S, p_S_progress) ## Random number of individuals become exposed based on FOI lambda
delta_S <- min(S, n_S_progress + rpois(seed_var1)) ## We make sure we don't remove all susceptibles


## Tricky bit -- deal with number becoming going through exposed and infected compartments
## EXPOSED
## This is the number of people who move from each stage
n_E_progress[] <- rbinom(E[i], p_E_progress)

aux_E[] <- (if (i == 1) delta_S  ## Number moving into each compartment is either number of lost S (first stage)
            else n_E_progress[i - 1]) - n_E_progress[i] ## Or the number moving into the class, less the number moving out
delta_E[] <- aux_E[i]

## INFECTED
## This is the number of people who move from each I stage
n_I_progress[] <- rbinom(I[i], p_I_progress) ## Number of new Infected this step

aux_I[] <- (if (i == 1) n_E_progress[k_E]  ## Number moving out of the last E class
            else n_I_progress[i - 1]) - n_I_progress[i] ## Or the number moving into the class, less the number moving out
delta_I[] <- aux_I[i]

## RECOVERED
delta_R <- n_I_progress[k_I]

## Now calculate the force of infection, which is the sum of contributions of all infected individuals modified by their stage-based modifier

#I_with_diff_trans[] <- trans_increase[i]*I[i]/N
#beta <- if(as.integer(step) >= tswitch2) beta3 else if(as.integer(step) >= tswitch1) beta2 else beta1
#beta <- max(betas[as.integer(step)], as.numeric(0))

lambda <- beta * sum(I[])/N


## Initial states are all zeroed as we will provide a state vector
## setting S and I based on the seeding model.
initial(S) <- S_ini
initial(E[]) <- 0
initial(E[1]) <- seed_value_var1
initial(I[]) <- 0
initial(R) <- 0
initial(cumu_inc) <- 0

## Parameters
#beta1 <- user(0.1)
#beta2 <- user(0.1)
#beta3 <- user(0.1)

#tswitch1 <- user(100)
#tswitch2 <- user(200)
#betas <- user()
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
#trans_increase[] <- user()
k_E <- user()
k_I <- user()
S_ini <- user()

## Dimensions
dim(E) <- k_E
dim(I) <- k_I

dim(n_E_progress) <- k_E
dim(aux_E) <- k_E
dim(delta_E) <- k_E

dim(n_I_progress) <- k_I
dim(aux_I) <- k_I
dim(delta_I) <- k_I


#dim(trans_increase) <- k_I
#dim(I_with_diff_trans) <- k_I