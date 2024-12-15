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

years <- seq(model$vec_date[1],2100,by=model$tstep)
temp <- T_RCP[years_MAGICC %in% years]

base_year <- model$vec_date[1]
T0 <- temp[1]
seq_i <- (years - base_year)/model$tstep

Dist4alpha_calibration <- function(par){
  Tinf  <- par[1]
  alpha <- par[2]
  traj_T <- Tinf - (Tinf - T0)*exp(-alpha * seq_i)
  distance <- (traj_T - temp)^2
  return(100000*sum(distance))
}

param_ini <- c(3,0.05)
for(i in 1:10){
  res_optim <- 
    optim(par = param_ini, # initial value of parameters
          fn = Dist4alpha_calibration, # function to optimize
          gr = NULL, # gradient of the function
          method = "BFGS", # optimization algorithm
          #method = "Nelder-Mead", # an alternative optimization algorithm
          control=list(trace=FALSE, # print progress
                       fnscale = 1, # fn scaling (-1: minimiz. -> optimiz.)
                       maxit=10000), # maximum number of iterations
          hessian = FALSE) # computes Hessian at the en dof optimiz.
  param_ini <- res_optim$par
}

Tinf  <- res_optim$par[1]
alpha <- res_optim$par[2]

traj_T <- Tinf - (Tinf - T0)*exp(-alpha * seq_i)


FILE = "/outputs/Figures/Figure_fit_alpha.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=6, height=3)

par(mfrow=c(1,1))
par(plt=c(.15,.95,.1,.95))

plot(years,temp,col="white",las=1,
     ylab="Temperature, in °C")
lines(years,traj_T,col="red",type="l",lwd=2)
lines(years,temp,col="black",lty=2,lwd=3)
lines(c(2020,2100),
      c(T_RCP[RCP_MAGICC$V1==2020],T_RCP[RCP_MAGICC$V1==2100]),col="blue",lty=3,lwd=2)

legend("topleft",
       legend=c("Average RCP45+RCP6",
                #expression(paste("Fitted parametric trajectory (",alpha,"=0.03, ",T[infinity],"=5.65)",sep="")),
                "Fitted parametric trajectory",
                "Linear trajectory"),
       lty=c(2,1,3),
       col=c("black","red","blue"),
       pch=c(NaN),
       lwd=c(3,2,2),
       seg.len = 3,
       bty = "n",cex=1)

grid()

dev.off()



