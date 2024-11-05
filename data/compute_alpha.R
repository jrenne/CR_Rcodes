# This script determines the value of the alpha parameter used in the 
# calibration exercise. This parameter determines the curvature of the
# conditional trajectory of atmospheric temperatures.

#RCP data
temp<-read.table("./data/mean_ssp.txt",header=TRUE)
temp <- temp[temp$year >= 2020,] # keep only dates after 2020

base_year <- 2020
T0 <- temp$mean[temp$year==base_year]
seq_i <- (temp$year - base_year)/model$tstep

Dist <- function(par){
  Tinf  <- par[1]
  alpha <- par[2]
  traj_T <- Tinf - (Tinf - T0)*exp(-alpha * seq_i)
  distance <- abs(traj_T - temp$mean)
  return(sum(distance))
}

res_optim <- 
  optim(par = c(4,0.04), # initial value of parameters
        fn = Dist, # function to optimize
        gr = NULL, # gradient of the function
        method = "BFGS", # optimization algorithm
        #method = "Nelder-Mead", # an alternative optimization algorithm
        control=list(trace=T, # print progress
                     fnscale = 1, # fn scaling (-1: minimiz. -> optimiz.)
                     maxit=100), # maximum number of iterations
        hessian = TRUE) # computes Hessian at the en dof optimiz.

Tinf  <- res_optim$par[1]
alpha <- res_optim$par[2]
traj_T <- Tinf - (Tinf - T0)*exp(-alpha * seq_i)

par(mfrow=c(1,1))
plot(temp$year,temp$mean,col="white",las=1)
lines(temp$year,traj_T,col="red",type="l",lwd=2)
lines(temp$year,temp$mean,col="black",lty=3,lwd=3)
grid()





