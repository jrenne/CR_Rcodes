

RCP_MAGICC <- read.csv("data/RCP_Mat_MAGICC.csv", header=FALSE)
RCP_ACE    <- read.csv("data/RCP_Mat_ACE.csv", header=FALSE)
# indic2020 <- which(RCP_MAGICC$V1==2020)
# RCP_MAGICC <- RCP_MAGICC[indic2020:dim(RCP_MAGICC)[1],]

H2500 <- model$horiz.2100 + (2500 - 2100)/model_sol$tstep
years <- seq(model_sol$vec_date[1],2300,by=model_sol$tstep)

#model_sol$parameters$m0 <- 1.924217
#model_sol$parameters$m0 <- 2.5

indicators <- which(RCP_MAGICC$V1 %in% years)

RCP30 <- RCP_MAGICC$V2[indicators]
RCP45 <- RCP_MAGICC$V3[indicators]
RCP60 <- RCP_MAGICC$V4[indicators]
RCP85 <- RCP_MAGICC$V5[indicators]

T_MAGICC_RCP30 <- RCP_MAGICC$V6[indicators]
T_MAGICC_RCP45 <- RCP_MAGICC$V7[indicators]
T_MAGICC_RCP60 <- RCP_MAGICC$V8[indicators]
T_MAGICC_RCP85 <- RCP_MAGICC$V9[indicators]

T_ACE_RCP30 <- RCP_ACE$V2
T_ACE_RCP45 <- RCP_ACE$V3
T_ACE_RCP60 <- RCP_ACE$V4
T_ACE_RCP85 <- RCP_ACE$V5

#Model-implied EV
EV<-EV.fct(model_sol,H2500)

Mat.trajectory <- EV$EX$M_at
Mat.trajectory <- RCP45[-1]
#Mat.trajectory <- RCP85[-1]

Mat.trajectory <- RCP85[-1]
Tat.MAGICC     <- T_MAGICC_RCP85
Tat.ACE        <- T_ACE_RCP85

Mat.trajectory <- RCP60[-1]
Tat.MAGICC     <- T_MAGICC_RCP60
Tat.ACE        <- T_ACE_RCP60

# Mat.trajectory <- RCP45[-1]
# Tat.MAGICC     <- T_MAGICC_RCP45
# Tat.ACE        <- T_ACE_RCP45

#model_sol_aux$tstep <- 5
res_CR    <- simul_TAT_condit_MAT(model_sol,Mat.trajectory)
res_CDICE <- simul_TAT_condit_MAT_CDICE(model_sol$tstep,Mat.trajectory,
                                        model_sol$vector.ini$ini_Tat,
                                        model_sol$vector.ini$ini_Tlo,
                                        model_sol$vector.ini$ini_F,
                                        model_sol$parameters)

Tat.nonlinear <- res_CR$Tat.nonlinear
Tat.linear    <- res_CR$Tat.linear

plot(years,Tat.nonlinear,type="l",lwd=2,
     xlab="year",ylab="Atm. temperature (Degrees Celsius)",
     ylim=c(1,1.2*max(Tat.nonlinear,Tat.MAGICC,Tat.ACE)),
     main=expression(paste("(c) Effect of linearization on atm. temperature trajectory")))
points(years,Tat.linear,pch=3,col="#00000044",lwd=2)
#lines(c(model_sol$vector.ini$ini_Tat,EV$EX$T_at),col="blue")

lines(years,Tat.MAGICC,col="red")
lines(RCP_ACE$V1,Tat.ACE,col="red",lty=2)

lines(years,res_CDICE$Tat,col="blue")

legend("bottomright",
       legend=c("Linearized","Non-linearized"),
       lty=c(NaN,1),
       col=c("#00000044","black"),
       pch=c(3,NaN),
       lwd=c(2,2),seg.len = 3,
       bty = "n",cex=1)


