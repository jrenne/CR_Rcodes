# ==============================================================================
# Figure with distributions of SLR
# ==============================================================================

# Prepare omega vectors (for pricing):
omega_ZCB <- matrix(0,model_sol$n.X)
a         <- matrix(0,length(model_sol$X),1)
HSL       <- which(model_sol$names.var.X=="H")
a[HSL]    <- 1

H <- model_sol$horiz.2100

# Load external data:
load("./data/sea_rcp.Rdat")
sea60.smooth <- predict(loess(sea.rcp[,3]~sea.rcp[,1]))
sea45.smooth <- predict(loess(sea.rcp[,2]~sea.rcp[,1]))

# For Fourier transform:
x <- exp(seq(-5,5,length.out = 1000)) #grid for Proposition 8 (Fourier)

values.of.sl <- seq(#round(model_sol$vector.ini$ini_H,1),
  0,
  round(EV$EX$H[H]+5*sqrt(EV$VX[[HSL]][H]),1),
  by=.01)

nb.values.variable <- 2000

cl <- makeCluster(number.of.cores)
registerDoParallel(cl)

save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))

clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))

# Compute P and Q proba:
Price.ZC <- varphi(model_sol,omega_ZCB,H = H)[[3]]
all.Probas.Q <- foreach(h = 1:H, .combine=cbind) %dopar% {
  SL.Q<-c(varphi.hat.fast(model_sol,omega = omega_ZCB,
                          H=h,x,a = a,
                          b=values.of.sl)/Price.ZC[h])
  SL.Q
}
all.Probas.P <- foreach(i = 1:H, .combine=cbind) %dopar% {
  SL.P <- fourier(model_sol,x,values.of.sl,i,HSL)
  SL.P
}

stopCluster(cl)
file.remove("outputs/toto.Rdata")


CI.P <- confidence_intervals_across_horizons(all.Probas.P,
                                             values.of.variable = values.of.sl,
                                             nb.values.variable = nb.values.variable,
                                             vector.of.CI = vector.of.CI)
CI.Q <- confidence_intervals_across_horizons(all.Probas.Q,
                                             values.of.variable = values.of.sl,
                                             nb.values.variable = nb.values.variable,
                                             vector.of.CI = vector.of.CI)
scale.sl.values <- CI.P$scale.variable.values
all.CI.P  <- CI.P$all.CI
all.CI.Q  <- CI.Q$all.CI
all.pdf.P <- CI.P$all.pdf
all.pdf.Q <- CI.Q$all.pdf
all.cdf.P <- CI.P$all.cdf
all.cdf.Q <- CI.Q$all.cdf


# Compute expected trajectories (under P and Q):
E.P    <- EV$EX$H[1:H]
E.Q    <- varphi.tilde(model_sol,a,H)[[1]]/
  varphi(model_sol,omega_ZCB,H)[[3]]

#Plots
FILE = paste("/outputs/Figures/Figure_SL_P_and_Q_vector_CI.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=7, height=5)

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE),
       widths=c(1,1), heights=c(1,1))
par(plt=c(.16,.96,.15,.85))

y.lim <- c(0.1,1.4)
x.lim <- c(2040,2100)

plot(model_sol$vec_date[2:(H+1)],E.P,
     xlim=x.lim,ylim=y.lim,
     cex.main=1.5,cex.axis=1.5,cex.lab=1.5,
     col="white",xlab="",ylab="in meters",las=1,
     main="(a) - Global sea level rise")
lines(model_sol$vec_date[2:(H+1)],
      E.P,lwd=2,col=P.col.line)
for(i in length(vector.of.CI):1){
  # P
  polygon(c(model_sol$vec_date[2:(H+1)],rev(model_sol$vec_date[2:(H+1)])),
          c(all.CI.P[1,,i],rev(all.CI.P[2,,i])),
          col=P.col,border = NA)
}

plot(model_sol$vec_date[2:(H+1)],E.Q,
     xlim=x.lim,ylim=y.lim,las=1,
     col="white",xlab="",ylab="in meters",
     cex.main=1.5,cex.axis=1.5,cex.lab=1.5,
     main="(b) - Global sea level rise (risk-adjusted)")
for(i in length(vector.of.CI):1){
  # Q
  polygon(c(model_sol$vec_date[2:(H+1)],rev(model_sol$vec_date[2:(H+1)])),
          c(all.CI.Q[1,,i],rev(all.CI.Q[2,,i])),
          col=Q.col,border = NA)
}
lines(model_sol$vec_date[2:(H+1)],
      E.Q,lwd=2,col=Q.col.line)
lines(model_sol$vec_date[2:(H+1)],
      E.P,lwd=2,col=P.col.line)

par(plt=c(.08,.98,.15,.85))

plot(scale.sl.values[2:(nb.values.variable+1)],all.pdf.P[,H],type="l",
     col=P.col.line,lwd=3,
     xlim=c(0,
            round(EV$EX$H[H]+3*sqrt(EV$VX[[HSL]][H]),1)),
     main="(c) - P.d.f. of global sea level rise in 2100",
     cex.main=1.5,cex.axis=1.5,cex.lab=1.5,
     xlab="",ylab="",yaxt="no")
lines(scale.sl.values[2:(nb.values.variable+1)],all.pdf.Q[,H],
      col=Q.col.line,lwd=3)
abline(v=scale.sl.values[which.min(abs(all.cdf.Q[,model_sol$horiz.2100]
                                       -0.5))],
       col=Q.col.line,lty=3,lwd=2)                                              #median of Q
abline(v=E.Q[H],col=Q.col.line,lty=1,lwd=1)                                    #mean of Q
abline(v=scale.sl.values[which.min(abs(all.cdf.P[,model_sol$horiz.2100]
                                       -0.5))],
       col=P.col.line,lty=3,lwd=2)                                              #median of P
abline(v=E.P[H],col=P.col.line,lty=1,lwd=1)                                    #mean of P

legend("topright",
       legend=c("Physical p.d.f.","Risk-adjusted p.d.f.","Mean","Median"),
       lty=c(1,1,1,3),
       col=c(P.col.line,Q.col.line,"black","black"),cex=1.5,
       lwd=c(2,2,1,2),bty = "n")
dev.off()
