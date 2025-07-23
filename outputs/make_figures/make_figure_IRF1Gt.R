# ==============================================================================
# FIGURE 2. Atmospheric temperature response to a carbon pulse
# Figure_IRF1GtC.pdf
# ==============================================================================

# CDICE: -----------------------------------------------------------------------
# Load RCPs (to have an emission scenario, as in Figure 2.2.3 of EPA, 2023):
RCP_MAGICC <- read.csv("data/RCP_Mat_MAGICC.csv", header=FALSE)
Emissions_CO2_RCP <- RCP_MAGICC[,11] # loads RCP4.5

for(dt in c(1,5)){
  MATeq <- 851
  MUOeq <- 628
  MLOeq <- 1323
  r1 <- MATeq/MUOeq
  r2 <- MUOeq/MLOeq
  b12 <- 0.054
  b11 <- - b12
  b23 <- 0.0082
  b13 <- 0
  b21 <- b12 * r1
  b22 <- - b21 - b23
  b31 <- 0
  b32 <- b23 * r2
  b33 <- - b32
  B <- matrix(c(b11,b12,b13,b21,b22,b23,b31,b32,b33),3,3)
  c1 <- .137
  c3 <- .73
  c4 <- .00689
  lambda <- 1.06
  F2XCO2 <- 3.45
  shock <- c(1,0,0)
  dates <- seq(2030,2300,by=dt)
  # Emissions:
  Emissions <- Emissions_CO2_RCP[which(RCP_MAGICC$V1 %in% dates)]
  if(dt > 1){
    for(iiii in 1:(dt-1)){
      Emissions <- Emissions +
        Emissions_CO2_RCP[which(RCP_MAGICC$V1 %in% dates) - iiii]
    }
  }
  all.M <- NULL
  all.F <- NULL
  all.TAT <- NULL
  all.TLO <- NULL
  all.M.shock <- NULL
  all.Forc.shock <- NULL
  all.TAT.shock <- NULL
  all.TLO.shock <- NULL
  t <- 0
  for(y in dates){
    t <- t + 1
    if(y==dates[1]){
      M <- matrix(c(MATeq,
                    MUOeq,MLOeq),ncol=1)
      M.shock <- matrix(c(MATeq+1,
                          MUOeq,MLOeq),ncol=1)
      TAT <- model_sol$vector.ini$ini_Tat
      TLO <- model_sol$vector.ini$ini_Tlo
      TAT.shock <- model_sol$vector.ini$ini_Tat
      TLO.shock <- model_sol$vector.ini$ini_Tlo
    }else{
      M <- (diag(3) + dt*B) %*% M + c(Emissions[t],0,0)
      M.shock <- (diag(3) + dt*B) %*% M.shock+ c(Emissions[t],0,0)
    }
    Forc <- F2XCO2 * log(M[1]/MATeq)/log(2)
    Forc.shock <- F2XCO2 * log(M.shock[1]/MATeq)/log(2)
    
    all.M         <- cbind(all.M,M)
    all.F         <- cbind(all.F,Forc)
    all.TAT       <- cbind(all.TAT,TAT)
    all.TLO       <- cbind(all.TLO,TLO)
    
    all.M.shock   <- cbind(all.M.shock,M.shock)
    all.Forc.shock   <- cbind(all.Forc.shock,Forc.shock)
    all.TAT.shock <- cbind(all.TAT.shock,TAT.shock)
    all.TLO.shock <- cbind(all.TLO.shock,TLO.shock)
    
    TAT_1 <- TAT # to be used un TLO
    TAT <- TAT + c1*dt * (Forc - lambda * TAT - c3 * (TAT - TLO))
    TLO <- TLO + c4*dt * (TAT_1 - TLO)
    
    TAT.shock_1 <- TAT.shock # to be used un TLO
    TAT.shock   <- TAT.shock + c1*dt * (Forc.shock - lambda * TAT.shock - 
                                          c3 * (TAT.shock - TLO.shock))
    TLO.shock <- TLO.shock + c4*dt * (TAT.shock_1 - TLO.shock)
  }
  if(dt==1){
    IRF_CDICE.dt1 <- all.TAT.shock - all.TAT
    IRF_CDICE.dt1_dates <- dates
  }
  if(dt==5){
    IRF_CDICE.dt5 <- all.TAT.shock - all.TAT
    IRF_CDICE.dt5_dates <- dates
  }
}

# ------------------------------------------------------------------------------

# CR: -----------------------------------------------------------------------
values_iniMat <- c(model_sol$vector.ini$ini_Mat,1800)
for(largeMat in c(TRUE,FALSE)){
  for(indic_lineariz in c(TRUE,FALSE)){
    for(dt in c(1,5)){
      MATeq <- ifelse(largeMat,values_iniMat[2],values_iniMat[1])
      MUOeq <- model_sol$vector.ini$ini_Mup
      MLOeq <- model_sol$vector.ini$ini_Mlo
      ksi1  <- model_sol$parameters$xi_1
      ksi2  <- model_sol$parameters$xi_2
      ksi3  <- model_sol$parameters$xi_3
      tau   <- model_sol$parameters$tau
      nu    <- model_sol$parameters$nu
      shock <- c(1,0,0)
      dates <- seq(2030,2300,by=dt)

      all.M          <- NULL
      all.F          <- NULL
      all.TAT        <- NULL
      all.TLO        <- NULL
      all.M.shock    <- NULL
      all.Forc.shock <- NULL
      all.TAT.shock  <- NULL
      all.TLO.shock  <- NULL
      ttt <- 0
      for(y in dates){
        ttt <- ttt + 1
        if(y==dates[1]){
          M <- matrix(c(MATeq,
                        MUOeq,MLOeq),ncol=1)
          M.shock <- matrix(c(MATeq+1,
                              MUOeq,MLOeq),ncol=1)
          TAT <- model_sol$vector.ini$ini_Tat
          TLO <- model_sol$vector.ini$ini_Tlo
          TAT.shock <- model_sol$vector.ini$ini_Tat
          TLO.shock <- model_sol$vector.ini$ini_Tlo
        }else{
          M <- (model_sol$varphi%^%dt) %*% M
          M.shock <- (model_sol$varphi%^%dt) %*% M.shock
        }
        if(indic_lineariz){
          Forc <- tau/model_sol$parameters$m0[ttt]/log(2)/model_sol$parameters$m_pi * 
            M[1]
          Forc.shock <- tau/model_sol$parameters$m0[ttt]/log(2)/model_sol$parameters$m_pi * 
            M.shock[1]
        }else{
          Forc <- tau * log(M[1]/MATeq)/log(2)
          Forc.shock <- tau * log(M.shock[1]/MATeq)/log(2)
        }
        
        all.M         <- cbind(all.M,M)
        all.F         <- cbind(all.F,Forc)
        all.TAT       <- cbind(all.TAT,TAT)
        all.TLO       <- cbind(all.TLO,TLO)
        
        all.M.shock   <- cbind(all.M.shock,M.shock)
        all.Forc.shock   <- cbind(all.Forc.shock,Forc.shock)
        all.TAT.shock <- cbind(all.TAT.shock,TAT.shock)
        all.TLO.shock <- cbind(all.TLO.shock,TLO.shock)
        
        TAT <- TAT + ksi1*dt * (Forc - tau/nu * TAT - ksi2 * (TAT - TLO))
        TLO <- TLO + ksi3*dt * (TAT - TLO)
        
        TAT.shock   <- TAT.shock + ksi1*dt * (Forc.shock - tau/nu * TAT.shock - 
                                                ksi2 * (TAT.shock - TLO.shock))
        TLO.shock <- TLO.shock + ksi3*dt * (TAT.shock - TLO.shock)
      }
      if((dt==1)&indic_lineariz&!largeMat){
        IRF_CR.linear.dt1 <- all.TAT.shock - all.TAT
        IRF_CR.linear.dt1_dates <- dates
      }
      if((dt==5)&indic_lineariz&!largeMat){
        IRF_CR.linear.dt5 <- all.TAT.shock - all.TAT
        IRF_CR.linear.dt5_dates <- dates
      }
      if((dt==1)&!indic_lineariz&!largeMat){
        IRF_CR.nonlinear.dt1 <- all.TAT.shock - all.TAT
        IRF_CR.nonlinear.dt1_dates <- dates
      }
      if((dt==5)&!indic_lineariz&!largeMat){
        IRF_CR.nonlinear.dt5 <- all.TAT.shock - all.TAT
        IRF_CR.nonlinear.dt5_dates <- dates
      }
      if((dt==1)&indic_lineariz&largeMat){
        IRF_CR.linear.largeMat.dt1 <- all.TAT.shock - all.TAT
        IRF_CR.linear.largeMat.dt1_dates <- dates
      }
      if((dt==5)&indic_lineariz&largeMat){
        IRF_CR.linear.largeMat.dt5 <- all.TAT.shock - all.TAT
        IRF_CR.linear.largeMat.dt5_dates <- dates
      }
      if((dt==1)&!indic_lineariz&largeMat){
        IRF_CR.nonlinear.largeMat.dt1 <- all.TAT.shock - all.TAT
        IRF_CR.nonlinear.largeMat.dt1_dates <- dates
      }
      if((dt==5)&!indic_lineariz&largeMat){
        IRF_CR.nonlinear.largeMat.dt5 <- all.TAT.shock - all.TAT
        IRF_CR.nonlinear.largeMat.dt5_dates <- dates
      }
    }
  }
}
# ------------------------------------------------------------------------------


# Load EPA data:
fair   <- read.csv("data/EPA_fair.csv")
hector <- read.csv("data/EPA_hector.csv")
magicc <- read.csv("data/EPA_magicc.csv")

lower.bound.01 <- fair$q01
upper.bound.99 <- fair$q99
lower.bound.05 <- fair$q05
upper.bound.95 <- fair$q95

shock <- 1 # Gt Carbon

indic.M_at <- which(model_sol$names.var.X=="M_at")
indic.Forc <- which(model_sol$names.var.X=="Forc")

model_noN       <- model_sol
model_noN$parameters$a_N <- 0
model_noN$parameters$b_N <- 0*model_noN$parameters$b_N
model_sol_noN <- model_solve(model_noN,indic_CRRA = FALSE)

model_sol_shock     <- model_sol
model_sol_noN_shock <- model_sol_noN

# implement shock:
model_sol_shock$vector.ini$ini_Mat <- model_sol_shock$vector.ini$ini_Mat + shock
model_sol_shock$vector.ini$ini_F   <- model_sol_shock$vector.ini$ini_F - 
  model_sol$A0.star.inf[5,6] * shock
model_sol_shock$X[indic.M_at] <- model_sol_shock$vector.ini$ini_Mat
model_sol_shock$X[indic.Forc] <- model_sol_shock$vector.ini$ini_F

model_sol_noN_shock$vector.ini$ini_Mat <- model_sol_shock$vector.ini$ini_Mat
model_sol_noN_shock$X[indic.M_at]      <- model_sol_shock$vector.ini$ini_Mat
model_sol_noN_shock$vector.ini$ini_F   <- model_sol_shock$vector.ini$ini_F
model_sol_noN_shock$X[indic.Forc]      <- model_sol_shock$vector.ini$ini_F

end.date <- 2300
h.end <- (end.date-model_sol$vec_date[1])/model_sol$tstep - 1

EV       <- EV.fct(model_sol,h=h.end)
EV_shock <- EV.fct(model_sol_shock,h=h.end)
IRF_TAT  <- EV_shock$EX$T_at - EV$EX$T_at

EV_noN       <- EV.fct(model_sol_noN,h=h.end)
EV_noN_shock <- EV.fct(model_sol_noN_shock,h=h.end)
IRF_noN_TAT  <- EV_noN_shock$EX$T_at - EV_noN$EX$T_at





# ------------------------------------------------------------------------------
# Plot ----
FILE = "/outputs/Figures/Figure_IRF1GtC.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=8, height=4)

par(mfrow=c(1,1))
par(plt=c(.2,1,.1,.95))

nf <- layout(
  matrix(c(1,2), ncol=2, byrow=TRUE), 
  widths=c(3,1.5), 
  heights=c(1)
)

plot(EV$date+5,c(0,IRF_TAT[1:(h.end-1)]),col="white",las=1,
     xlim=c(2025,2300),
     ylim=c(0,1.5*max(upper.bound.99)),
     xlab="",
     ylab="",main="")
grid()

title(ylab='Temperature Anomaly from 1GtC, in °C', line=5,
      cex.lab=1)

polygon(x=c(fair$year,rev(fair$year)),c(lower.bound.01,rev(upper.bound.99)),
        col='grey90',border=NaN)
polygon(x=c(fair$year,rev(fair$year)),c(lower.bound.05,rev(upper.bound.95)),
        col='grey80',border=NaN)

lines(EV$date,c(0,IRF_TAT[1:(h.end-1)]),lwd=1,pch=3,type="b")
lines(EV$date,c(0,IRF_noN_TAT[1:(h.end-1)]),lwd=1)
lines(fair$year,fair$mean,lwd=2,col="dark grey")
lines(hector$year,hector$temp.delta,col="#E69F00",lty=3,lwd=2)
lines(magicc$year,magicc$temp.delta,col="#56B4E9",lty=2,lwd=2)
lines(IRF_CDICE.dt5_dates,IRF_CDICE.dt5,col="red",lwd=2,lty=4)

# Load Traeger (2023) responses:
IRF_Traeger_5y <- read.csv("data/IRF_Traeger_5y.csv", header=FALSE)
n <- length(IRF_Traeger_5y$V1)
dates <- seq(2030,2030+5*(n-1),by=5)
lines(dates,IRF_Traeger_5y$V1/100,col="#006400",lty=4)
lines(dates,IRF_Traeger_5y$V2/100,col="#006400",lty=5)

plot.new()

par(plt=c(.1,.95,.25,.85))

legend("topleft",
       legend=c("CR (with N)","CR (w/o N)","HECTOR 2.5","MAGICC 7.5.3",
                "FaIR 1.6.2 (mean)","CDICE (5yr time step)",
                "ACE-DICE (5yr time step)","ACE-Joos (5yr time step)"),
       col=c("black","black","#E69F00", "#56B4E9","dark grey","red",
             "#006400","#006400"),
       lty=c(NaN,1,3,2,1,4,4,5),
       lwd=c(1,1,2,2,2,2,1,1),
       pch=c(3,NaN,NaN,NaN,NaN,NaN,NaN,NaN),
       bty = "n",cex=1,
       bg="white",
       seg.len = 3)
dev.off()
