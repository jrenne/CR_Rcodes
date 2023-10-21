#radiative forcings w/o approx
#2040
H40<-model_sol$horiz.2100-12

sigma_F <- .1

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

Mat.distr.2040 <-fourier(model_sol,x,gridMat,H40,6)
Mat.pdf.2040   <- diff(Mat.distr.2040)
Mat.distr.2100 <-fourier(model_sol,x,gridMat,H,6)
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
pdf(file=paste(getwd(),FILE,sep=""),pointsize=8,width=6, height=6)
par(mfrow=c(2,1))
par(plt=c(.15,.95,.2,.85))
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
       lwd=c(2,2),
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
       legend=c("2040","2100"),
       lty=c(1,1),
       col=c(P.col.line,P.col.line),
       lwd=c(3,1),
       bty = "n",cex=1)

dev.off()