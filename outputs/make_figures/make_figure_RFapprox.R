# ==============================================================================
# Figure illustrating the linearization of radiative forcings' equation
# ==============================================================================

interm.year <- 2050
H40 <- which(model_sol$vec_date==interm.year)-1

sigma_F <- .001

# For Fourier-transform computation:
x <- exp(seq(-10,10,length.out = 2000)) # grid for Proposition 8 (Fourier)

param <- model$parameters

gridMat<-seq(900,1800,by=1)

# Var_eta2040<-EV$VX[[model_sol$n.Z+2]][H40]*
#   ((1-param$Phi[2,2]^2)^(1/2)*param$sigma_F)^2
# Var_eta2100<-EV$VX[[model_sol$n.Z+2]][H]*
#   ((1-param$Phi[2,2]^2)^(1/2)*param$sigma_F)^2
Var_eta2040 <- sigma_F^2
Var_eta2100 <- sigma_F^2

nlinF <- matrix(0,length(gridMat),1)
nlinF[1:length(nlinF)] <- model_sol$parameters$tau/log(2)*
  log(gridMat[1:length(gridMat)]/model_sol$parameters$m_pi)

linF  <-matrix(0,length(gridMat),1)
linF[1:length(linF)]  <- model_sol$parameters$tau/log(2)*
  (log(model_sol$parameters$m0)+
     (gridMat[1:length(gridMat)]/model_sol$parameters$m_pi-
        model_sol$parameters$m0)/(model_sol$parameters$m0))

indic.Mat <- which(model_sol$names.var.X=="M_at")

Mat.distr.2040 <-fourier(model_sol,x,gridMat,H40,indic.Mat)
Mat.pdf.2040   <- diff(Mat.distr.2040)
Mat.distr.2100 <-fourier(model_sol,x,gridMat,model_sol$horiz.2100,indic.Mat)
Mat.pdf.2100   <- diff(Mat.distr.2100)

#x=M_at, y=F
gridF    <- seq(0.00,8,by=.002)
gridMat_1<- gridMat[-1]
nF       <- length(gridF)
nMat     <- length(gridMat_1)
mx.F     <- t(matrix(gridF,nF,nMat))
mx.Mat   <-   matrix(gridMat_1,nMat,nF)

#2040
pdf.cond <- 1/sqrt(2*pi*Var_eta2040)*
  exp(-(mx.F-model_sol$parameters$tau/log(2)*
          (log(model_sol$parameters$m0)+
             (mx.Mat/model_sol$parameters$m_pi-model_sol$parameters$m0)/
             (model_sol$parameters$m0)))^2/(2*Var_eta2040)
  )
pdf_xy <-  matrix(Mat.pdf.2040,ncol=nF,nrow = nMat)*pdf.cond
p <- .9
res <- make_confidence_area(gridMat_1,gridF,pdf_xy,p)
ca.Mat   <- res$x.polygon
ca.F     <- res$y.polygon


#2100
pdf.cond.2100 <- 1/sqrt(2*pi*Var_eta2100)*
  exp(-(mx.F-model_sol$parameters$tau/log(2)*
          (log(model_sol$parameters$m0)+
             (mx.Mat/model_sol$parameters$m_pi-model_sol$parameters$m0)/
             (model_sol$parameters$m0)))^2/(2*Var_eta2100)
  )
pdf_xy.2100 <-  matrix(Mat.pdf.2100,ncol=nF,nrow = nMat)*pdf.cond.2100
p <- .9
res.2100 <- make_confidence_area(gridMat_1,gridF,pdf_xy.2100,p)
ca.Mat.2100   <- res.2100$x.polygon
ca.F.2100     <- res.2100$y.polygon


#Plot
FILE = paste("/outputs/Figures/Figure_F_approx.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=6, height=6)
par(mfrow=c(3,1))
par(plt=c(.1,.95,.25,.85))
plot(ca.Mat, ca.F,col="white",ylim=c(1,6),xlim=range(gridMat),type="l",
     las=1,xlab=expression(paste("M"[AT]," (GtC)")),
     ylab=expression(paste("FCO"[2]," (in Wm-2)")),
     main=expression(paste("(a) Relationship between radiative forcings and atmospheric carbon concentration")))
# polygon(ca.Mat,ca.F,
#         col=adjustcolor("grey17", alpha.f = 0.15), border = NaN)
# polygon(ca.Mat.2100,ca.F.2100,
#         col=adjustcolor("grey17", alpha.f = 0.3), border = NaN)
lines(gridMat,nlinF,
      col="black",lwd=2)
lines(gridMat,linF,
      col="black",lwd=2,lty=2)
legend("topleft",
       legend=c("Linearized","Non-linearized"),
       lty=c(2,1),
       col=c("black","black"),
       lwd=c(2,2),seg.len = 3,
       bty = "n",cex=1)
# legend("bottomright",
#        title = "95% confidence interval:",
#        legend=c("2040","2100"),
#        fill=c(adjustcolor("grey17", alpha.f = 0.15), 
#               adjustcolor("grey17", alpha.f = 0.3)),
#        density=c(NaN, NaN),
#        bty = "n",cex=1)

plot(gridMat[-1],Mat.pdf.2040,type="l",yaxt="n",
     las=1,ylab="",xlab=expression(paste("M"[AT]," (GtC)")),
     col=P.col.line,lwd=3,
     main=expression(paste("(b) Atmospheric carbon concentration p.d.f.")))
lines(gridMat[-1],Mat.pdf.2100,
      col=P.col.line,lwd=1,lty=1)
legend("topleft",
       legend=c(interm.year,"2100"),
       lty=c(1,1),
       col=c(P.col.line,P.col.line),
       lwd=c(3,1),seg.len = 3,
       bty = "n",cex=1)

# Panel (c) illustrate effects on temperature trajectory -----------------------

H2500 <- model$horiz.2100 + (2500 - 2100)/model_sol$tstep
#Model-implied EV
EV<-EV.fct(model_sol,H2500)
Mat.trajectory <- EV$EX$M_at

param <- model_sol$parameters
H2100                <- model_sol$horiz.2100
f_ex                 <- matrix(rep(param$phi_0,H2500),H2500,1)
f_ex[1:H2100]        <- f_ex[1:H2100] +
  (1/H2100)*(param$phi_1-param$phi_0)*((1:H2100)-1)
f_ex[(H2100+1):H2500] <- f_ex[(H2100+1):H2500] + (param$phi_1-param$phi_0)


Mat <- model_sol$vector.ini$ini_Mat
Mup <- model_sol$vector.ini$ini_Mup
Mlo <- model_sol$vector.ini$ini_Mlo
Tat_1 <- model_sol$vector.ini$ini_Tat
Tlo_1 <- model_sol$vector.ini$ini_Tlo
F_1   <- model_sol$vector.ini$ini_F

xi_1 <- model_sol$parameters$xi_1
xi_2 <- model_sol$parameters$xi_2
xi_3 <- model_sol$parameters$xi_3

tau <- model_sol$parameters$tau
nu  <- model_sol$parameters$nu

m_pi <- model_sol$parameters$m_pi
m0   <- model_sol$parameters$m0

phi_0 <- model_sol$parameters$phi_0
phi_1 <- model_sol$parameters$phi_1

# phi11 <- model_sol$parameters$varphi_11
# phi12 <- model_sol$parameters$varphi_12
# phi13 <- model_sol$parameters$varphi_13
# phi21 <- model_sol$parameters$varphi_21
# phi22 <- model_sol$parameters$varphi_22
# phi23 <- model_sol$parameters$varphi_23
# phi31 <- model_sol$parameters$varphi_31
# phi32 <- model_sol$parameters$varphi_32
# phi33 <- model_sol$parameters$varphi_33
# Phi <- matrix(c(phi11,phi21,phi31,
#                 phi12,phi22,phi32,
#                 phi13,phi23,phi33),3,3)

## check:
#1 - phi12 * model_sol$parameters$mateq/model_sol$parameters$mueq - phi23
#phi22

for(linearized in 0:1){
  Tat_1 <- model_sol$vector.ini$ini_Tat
  Tlo_1 <- model_sol$vector.ini$ini_Tlo
  F_1   <- model_sol$vector.ini$ini_F
  all.Tat <- Tat_1
  all.F   <- F_1
  for(t in 1:length(Mat.trajectory)){
    Tat <- Tat_1 + xi_1 * (F_1 - tau/nu * Tat_1 - xi_2 * (Tat_1 - Tlo_1))
    Tlo <- Tlo_1 + xi_3 * (Tat_1 - Tlo_1)

    if(linearized==0){
      FF <- tau * log(Mat.trajectory[t]/m_pi)/log(2) + f_ex[t]
    }else{
      FF <- tau * log(m0)/log(2) +
        tau/(log(2)*m0)*(Mat.trajectory[t]/m_pi - m0) + f_ex[t]
    }

    all.Tat <- c(all.Tat,Tat)
    all.F   <- c(all.F,FF)

    Tat_1 <- Tat
    Tlo_1 <- Tlo
    F_1   <- FF
  }
  if(linearized==1){
    all.Tat.linear <- all.Tat
    all.F.linear   <- all.F
  }else{
    all.Tat.nonlinear <- all.Tat
    all.F.nonlinear   <- all.F
  }
}
# 
# years <- seq(model_sol$vec_date[1],2500,by=model_sol$tstep)
# plot(years,all.Tat.nonlinear,type="l",lwd=2,
#      xlab="year",ylab="Atm. temperature (Degrees Celsius)",
#      ylim=c(1,3.5),
#      main=expression(paste("(c) Effect of linearization on atm. temperature trajectory")))
# points(years,all.Tat.linear,pch=3,col="#00000044",lwd=2)
# #lines(c(model_sol$vector.ini$ini_Tat,EV$EX$T_at),col="blue")
# 
# legend("bottomright",
#        legend=c("Linearized","Non-linearized"),
#        lty=c(NaN,1),
#        col=c("#00000044","black"),
#        pch=c(3,NaN),
#        lwd=c(2,2),seg.len = 3,
#        bty = "n",cex=1)
# 
# # plot(all.F.nonlinear,type="l")
# # lines(all.F.linear,col="red")
# # lines(c(model_sol$vector.ini$ini_F,EV$EX$Forc),col="blue")

dev.off()

