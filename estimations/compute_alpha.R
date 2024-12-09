# This script determines the value of the alpha parameter used in the 
# calibration exercise. This parameter determines the curvature of the
# conditional trajectory of atmospheric temperatures.

#RCP data
temp<-read.table("./data/mean_ssp.txt",header=TRUE)
temp <- temp[temp$year >= 2020,] # keep only dates after 2020

RCP_MAGICC <- read.csv("data/RCP_Mat_MAGICC.csv", header=FALSE)
years_MAGICC <- RCP_MAGICC$V1
T_RCP_45 <- RCP_MAGICC$V7
T_RCP_60 <- RCP_MAGICC$V8
T_RCP <- .5*(T_RCP_45 + T_RCP_60)

years <- seq(model_sol$vec_date[1],2100,by=model_sol$tstep)
temp <- T_RCP[years_MAGICC %in% years]

base_year <- model_sol$vec_date[1]
T0 <- temp[1]
seq_i <- (years - base_year)/model_sol$tstep

Dist4alpha_calibration <- function(par){
  Tinf  <- par[1]
  alpha <- par[2]
  traj_T <- Tinf - (Tinf - T0)*exp(-alpha * seq_i)
  distance <- abs(traj_T - temp)
  return(sum(distance))
}

res_optim <- 
  optim(par = c(8,0.04), # initial value of parameters
        fn = Dist4alpha_calibration, # function to optimize
        gr = NULL, # gradient of the function
        method = "BFGS", # optimization algorithm
        #method = "Nelder-Mead", # an alternative optimization algorithm
        control=list(trace=T, # print progress
                     fnscale = 1, # fn scaling (-1: minimiz. -> optimiz.)
                     maxit=1000), # maximum number of iterations
        hessian = TRUE) # computes Hessian at the en dof optimiz.

Tinf  <- res_optim$par[1]
alpha <- res_optim$par[2]
traj_T <- Tinf - (Tinf - T0)*exp(-alpha * seq_i)

par(mfrow=c(1,1))
plot(years,temp,col="white",las=1)
lines(years,traj_T,col="red",type="l",lwd=2)
lines(years,temp,col="black",lty=3,lwd=3)
grid()





