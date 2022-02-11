gen_sir <- odin.dust::odin_dust("../code/sirtest.R")
k_E <- 2
k_I <- 2
S_ini = 1000000
dt <- 0.25
beta_step <- c(rep(0.25,75/dt),rep(0.1,50/dt),rep(0.2,200/dt))
beta_step <- rep(0.25, 100000)
sir_model <- gen_sir$new(pars = list(dt = dt,
                                     S_ini = S_ini,
                                     beta_step = beta_step,
                                     gamma = 1/7,
                                     sigma=1/3,
                                     trans_increase=rep(1, k_I),
                                     k_E=k_E,
                                     k_I=k_I,
                                    # seed_start_var1=0,
                                    # seed_duration_var1=0.01,
                                     seed_value_var1=100),
                         step = 1,
                         n_particles = 10L,
                         n_threads = 4L,
                         seed = 1L)

n_steps <- 200/dt
n_particles <- 10L
x <- array(NA, dim = c(sir_model$info()$len, n_particles, n_steps))
for (t in seq_len(n_steps)) {
    x[ , , t] <- sir_model$run(t)
}
time <- x[1, 1, ]
x <- x[-1, , ]
par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
cols <- c(S = "#8c8cd9", E="black", I = "#cc0044", R = "#999966")
matplot(time, t(x[1, , ]), type = "l",
        xlab = "Time", ylab = "Number of individuals",
        col = cols[["S"]], lty = 1, ylim = range(x/S_ini))
matlines(time, t(x[2, , ]/S_ini), col = cols[["R"]], lty = 1)
E <- x[4,,]
if(k_E > 1){
    for(i in 5:(4+k_E-1)){
     E <- E + x[i,,]   
    }
}
matlines(time, t(E/S_ini), col = cols[["E"]], lty = 1)
I <- x[4+k_E,,]
if(k_I > 1){
    for(i in (4+k_E):(4+k_E+k_I-1)){
        I <- I + x[i,,]   
    }
}
matlines(time, t(I/S_ini), col = cols[["I"]], lty = 1)

matlines(time, t(x[3,,]/S_ini),col="orange",lty=1)

legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")

print(mean(x[2,,dim(x)[3]]))

#matplot(time, t(x[3,,]/S_ini),lty=1,type = "l",
#        xlab = "Time", ylab = "Number of individuals",ylim=range(x[3,,])/S_ini)

matplot(seq(1,200,by=1), apply(x[3,,]/S_ini, 1, function(y) rowSums(matrix(y, ncol=1/dt,byrow=TRUE))),lty=1,type = "l",
        xlab = "Time", ylab = "Number of individuals",ylim=c(0,0.02))
#lines(diff(c(0,sol1[,ncol(sol1)]/sum(N0_1))),type='l',col="orange",lwd=5)

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
#par(mfrow=c(2,1))
#plot(dgamma(0:25, shape=k_E,rate=k_E/3),type='l')
#plot(dgamma(0:25, shape=k_I,rate=k_I/7),type='l')
#par(mfrow=c(1,1))
#mean(rgamma(10000, shape=k_I,rate=k_I/7))
#pgamma(0,shape=k_I,scale=k_I*(1/7),lower.tail=FALSE)

