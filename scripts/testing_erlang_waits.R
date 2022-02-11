general_sir_erlang <- function(t,y, pars){
    beta <- pars[1] ## Overall transmission rate
    Ta <- pars[2] ## Incubation period
    Tg <- pars[3] ## Infectious period
    
    NE <- pars[4]
    NI <- pars[5]
    
    N <- sum(y)
    
    S <- y[1]
    Es <- y[2:(2+NE-1)]
    Is <- y[(2+NE):(2+NE+NI-1)]
    R <- y[length(y)]
    
    dEs <- numeric(NE)
    dIs <- numeric(NI)
    
    dS <- -beta*S*sum(Is)/N
    dEs[1] <- beta*S*sum(Is)/N - NE*Ta*Es[1]
    
    if(NE >= 2){
        for(i in 2:NE){
            dEs[i] <- NE*Ta*Es[i-1] - NE*Ta*Es[i]
        }
    }
    
    dIs[1] <- NE*Ta*Es[length(Es)] - NI*Tg*Is[1]
    
    if(NI >= 2){
        for(i in 2:NI){
            dIs[i] <- NI*Tg*Is[i-1] - NI*Tg*Is[i]
        }
    }
    dR <- NI*Tg*Is[length(Is)]
    
    return(list(c(dS,dEs,dIs,dR)))
}

general_sir <- function(t,y, pars){
    beta <- pars[1] ## Overall transmission rate
    Ta <- pars[2] ## Incubation period
    Tg <- pars[3] ## Infectious period
    
    N <- sum(y)
    
    S <- y[1]
    E <- y[2]
    I <- y[3]
    R <- y[4]
    
    dS <- -beta*S*I/N
    dE <- beta*S*I/N - Ta*E
    dI <- Ta*E - Tg*I
    dR <- Tg*I
    
    return(list(c(dS,dE,dI,dR)))
}

NE <- 2
NI <- 2
pars <- c(0.25,1/3,1/7,NE,NI)

Es <- rep(0, NE)
Es[1] <- 100

Is <- rep(0, NI)

N0_1 <- c(1000000,Es, Is, 0)
N0_2 <- c(1000000,1000,0,0)
times <- seq(1,200,by=1)
sol1 <- deSolve::ode(y = N0_1,t = times,parms=pars,func=general_sir_erlang)
sol2 <- deSolve::ode(y = N0_2,t = times,parms=pars,func=general_sir)

#plot(sol1[,ncol(sol1)]/sum(N0_1),type='l',col="blue",ylim=c(0,1))
#lines(sol2[,ncol(sol2)]/sum(N0_2),col="red")


#plot(diff(c(0,sol1[,ncol(sol1)]/sum(N0_1))),type='l',col="blue",ylim=c(0,0.02))
#lines(diff(c(0,sol2[,ncol(sol2)]/sum(N0_2))),col="red")
