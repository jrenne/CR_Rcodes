
# Prepare omega vectors (for pricing):
omega_ZCB <- matrix(0,model_sol$n.X)
omega_T.at <- omega_ZCB
omega_T.at[which(model_sol$names.var.X=="T_at")] <- 1

H <- model_sol$horiz.2100

load("./data/mat.Rdat")
par(plt=c(.1,.9,.1,.9))

x <- exp(seq(-10,10,length.out = 1000)) #grid for Proposition 8 (Fourier)

MAT<-which(model$names.var.X=="M_at")

values.of.emissions <- seq(round(model_sol$vector.ini$ini_Mat)-100,
                           round(EV$EX$M_at[H]+5*sqrt(EV$VX[[MAT]][H]))+100,
                           by=20)

a   <- matrix(0,length(model_sol$X),1)
a[MAT]<-1

#if(rerun==1){
cl <- makeCluster(number.of.cores)
registerDoParallel(cl)

save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))

clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))

# Compute P and Q proba:
Price.ZC <- varphi(model_sol,omega_ZCB,H = H)[[3]]
all.Probas.Q <- foreach(h = 1:H, .combine=cbind) %dopar% {
  Emissions.Q<-c(varphi.hat.fast(model_sol,omega = omega_ZCB,
                                 H=h,x,a = a,
                                 b=values.of.emissions)/Price.ZC[h])
  Emissions.Q
}
all.Probas.P <- foreach(i = 1:H, .combine=cbind) %dopar% {
  Emissions.P <- fourier(model_sol,x,values.of.emissions,i,MAT)
  Emissions.P
}

stopCluster(cl)

CI.P <- confidence_intervals_across_horizons(all.Probas.P,
                                             values.of.variable = values.of.emissions,
                                             nb.values.variable = nb.values.variable,
                                             vector.of.CI = vector.of.CI)
CI.Q <- confidence_intervals_across_horizons(all.Probas.Q,
                                             values.of.variable = values.of.emissions,
                                             nb.values.variable = nb.values.variable,
                                             vector.of.CI = vector.of.CI)
scale.emissions.values <- CI.P$scale.variable.values
all.CI.P  <- CI.P$all.CI
all.CI.Q  <- CI.Q$all.CI
all.pdf.P <- CI.P$all.pdf
all.pdf.Q <- CI.Q$all.pdf
all.cdf.P <- CI.P$all.cdf
all.cdf.Q <- CI.Q$all.cdf


# Compute expected trajectories (under P and Q):
E.P     <- EV$EX$M_at[1:H]
E.Q     <- varphi.tilde(model_sol,a,H)[[1]]/varphi(model_sol,omega_ZCB,H)[[3]]

#Plot
FILE = paste("/outputs/Figures/Figure_Mat_P_and_Q_vector_CI.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=9, height=6)

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE),
       widths=c(1,1,2), heights=c(1,1,2))
par(plt=c(.15,.95,.15,.85))

y.lim <- c(min(values.of.emissions),max(values.of.emissions))
x.lim <- c(2040,2100)

plot(model_sol$vec_date[2:(H+1)],E.P,
     xlim=x.lim,ylim=y.lim,
     cex.main=1.5,cex.axis=1.5,cex.lab=1.5,
     col="white",xlab="",ylab="",las=1,
     main="(a) - Carbon concentration in the atm.")

for(i in length(vector.of.CI):1){
  # P
  polygon(c(model_sol$vec_date[2:(H+1)],
            rev(model_sol$vec_date[2:(H+1)])),
          c(all.CI.P[1,1:H,i],rev(all.CI.P[2,1:H,i])),
          col=P.col,border = NA)
}
lines(model_sol$vec_date[2:(H+1)],
      E.P,lwd=2,col=P.col.line)
lines(mat[,1],mat[,2],lty=2,lwd=2,col="lightsteelblue4")
lines(mat[,1],mat[,3],lty=2,lwd=2,col="grey28")
legend("topright",
       legend=c("RCP4.5","RCP6.0"),
       lty=c(2,2),
       col=c("lightsteelblue4","grey28"),
       cex=1.5,
       lwd=c(2,2),bty = "n")

plot(model_sol$vec_date[2:(H+1)],E.Q,
     xlim=x.lim,ylim=y.lim,las=1,
     col="white",xlab="",ylab="",
     cex.main=1.5,cex.axis=1.5,cex.lab=1.5,
     main="(b) - Carbon concentration in the atm. (risk-adjusted)")
for(i in length(vector.of.CI):1){
  # Q
  polygon(c(model_sol$vec_date[2:(H+1)],
            rev(model_sol$vec_date[2:(H+1)])),
          c(all.CI.Q[1,1:H,i],rev(all.CI.Q[2,1:H,i])),
          col=Q.col,border = NA)
}
lines(model_sol$vec_date[2:(H+1)],
      E.Q,lwd=2,col=Q.col.line)
lines(model_sol$vec_date[2:(H+1)],
      E.P,lwd=2,col=P.col.line)
plot(scale.emissions.values[2:(nb.values.variable+1)],all.pdf.P[,H],type="l",
     col=P.col.line,lwd=3,
     xlim=c(EV$EX$M_at[H]-4*sqrt(EV$VX[[MAT]][H]),
            EV$EX$M_at[H]+4*sqrt(EV$VX[[MAT]][H])),
     main="(c) - P.d.f. of emissions in 2100",
     cex.main=1.5,cex.axis=1.5,cex.lab=1.5,
     xlab="",ylab="",yaxt="no")
lines(scale.emissions.values[2:(nb.values.variable+1)],all.pdf.Q[,H],
      col=Q.col.line,lwd=3)
abline(v=scale.emissions.values[which.min(abs(all.cdf.Q[,model_sol$horiz.2100]
                                              -0.5))],
       col=Q.col.line,lty=3,lwd=2)                                              #median of Q
abline(v=E.Q[H],col=Q.col.line,lty=1,lwd=1)                                    #mean of Q
abline(v=scale.emissions.values[which.min(abs(all.cdf.P[,model_sol$horiz.2100]
                                              -0.5))],
       col=P.col.line,lty=3,lwd=2)                                              #median of P
abline(v=E.P[H],col=P.col.line,lty=1,lwd=1)                                    #mean of P

legend("topright",
       legend=c("Physical p.d.f.","Risk-adjusted p.d.f.","Mean","Median"),
       lty=c(1,1,1,3),
       col=c(P.col.line,Q.col.line,"black","black"),cex=1.5,
       lwd=c(2,2,1,2),bty = "n")
dev.off()
