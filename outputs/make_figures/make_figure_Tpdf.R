# ==============================================================================
# Figure with temperature distributions
# ==============================================================================

H <- model_sol$horiz.2100

# Prepare omega vectors (for pricing):
omega_ZCB <- matrix(0,model_sol$n.X)
omega_T.at <- omega_ZCB
omega_T.at[which(model_sol$names.var.X=="T_at")] <- 1
TAT <- which(model_sol$names.var.X=="T_at")

#RCP data
temp<-read.table("./data/mean_ssp.txt",header=TRUE)
temp_graph<-temp[3:11,]

x <- exp(seq(-5,5,length.out = 1000)) #grid for Proposition 8 (Fourier)

#distribution and expected temperature
values.of.temperatures <- seq(round(model_sol$vector.ini$ini_Tat),
                              round(EV$EX$T_at[H]+5*sqrt(EV$VX[[TAT]][H])),
                              by=.1)

cl <- makeCluster(number.of.cores)
registerDoParallel(cl)

save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))

clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))

# Compute P and Q proba:
Price.ZC <- varphi(model_sol,omega_ZCB,H = H)[[3]]
all.Probas.Q <- foreach(h = 1:H, .combine=cbind) %dopar% {
  Temperature.Q<-c(varphi.hat.fast(model_sol,omega = omega_ZCB,
                                   H=h,x,a = omega_T.at,
                                   b=values.of.temperatures)/Price.ZC[h])
  Temperature.Q
}
all.Probas.P <- foreach(h = 1:H, .combine=cbind) %dopar% {
  Temperature.P <- fourier(model_sol,x,values.of.temperatures,h,TAT)
  Temperature.P
}

stopCluster(cl)
file.remove("outputs/toto.Rdata")


# Compute confidence intervals:
CI.P <- confidence_intervals_across_horizons(all.Probas.P,
                                                 values.of.variable = values.of.temperatures,
                                                 nb.values.variable = nb.values.variable,
                                                 vector.of.CI = vector.of.CI)
CI.Q <- confidence_intervals_across_horizons(all.Probas.Q,
                                                 values.of.variable = values.of.temperatures,
                                                 nb.values.variable = nb.values.variable,
                                                 vector.of.CI = vector.of.CI)
scale.temperatures.values <- CI.P$scale.variable.values
all.CI.P  <- CI.P$all.CI
all.CI.Q  <- CI.Q$all.CI
all.pdf.P <- CI.P$all.pdf
all.pdf.Q <- CI.Q$all.pdf
all.cdf.P <- CI.P$all.cdf
all.cdf.Q <- CI.Q$all.cdf

# Compute expected trajectories (under P and Q):
EV      <- EV.fct(model_sol,h=H)
ET.P    <- EV$EX$T_at[1:H]
ET.Q    <- varphi.tilde(model_sol,omega_T.at,H)[[1]]/
  varphi(model_sol,omega_ZCB,H)[[3]]

#Plots
FILE = paste("/outputs/Figures/Figure_Tat_P_and_Q_vector_CI.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=9, height=6)

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE),
       widths=c(1,1,2), heights=c(1,1,2))
par(plt=c(.15,.95,.15,.85))

y.lim <- c(.5,6)
x.lim <- c(2035,2100)

plot(model_sol$vec_date[2:(H+1)],ET.P,
     xlim=x.lim,ylim=y.lim, cex.main=1.5,cex.axis=1.5,cex.lab=1.5,
     col="white",xlab="",ylab="Degrees Celsius",las=1,
     main="(a) - Atmosph. temperature")

for(i in length(vector.of.CI):1){
  # P
  polygon(c(model_sol$vec_date[2:(H+1)],rev(model_sol$vec_date[2:(H+1)])),
          c(all.CI.P[1,,i],rev(all.CI.P[2,,i])),
          col=P.col,border = NA)
}
lines(model_sol$vec_date[2:(H+1)],
      ET.P,lwd=2,col=P.col.line)
lines(temp_graph[,1],temp_graph[,5],lty=2,lwd=1,col="lightsteelblue4")
lines(temp_graph[,1],temp_graph[,4],lty=2,lwd=1,col="lightsteelblue4")
lines(temp_graph[,1],temp_graph[,2],lty=2,lwd=2,col="grey28")
legend("topright",
       legend=c("RCP4.5+RCP6.0","+/- 2 std"),
       lty=c(2,2),
       col=c("grey28","lightsteelblue4"),cex=1.5,
       lwd=c(2,1),bty = "n")

plot(model_sol$vec_date[2:(H+1)],ET.Q,
     xlim=x.lim,ylim=y.lim,las=1,
     col="white",xlab="",ylab="Degrees Celsius", cex.main=1.5,cex.axis=1.5,
     cex.lab=1.5,
     main="(b) - Atmosph. temperature (risk-adjusted)")
for(i in length(vector.of.CI):1){
  # Q
  polygon(c(model_sol$vec_date[2:(H+1)],rev(model_sol$vec_date[2:(H+1)])),
          c(all.CI.Q[1,,i],rev(all.CI.Q[2,,i])),
          col=Q.col,border = NA)
}
lines(model_sol$vec_date[2:(H+1)],
      ET.Q,lwd=2,col=Q.col.line)
lines(model_sol$vec_date[2:(H+1)],
      ET.P,lwd=2,col=P.col.line)

plot(scale.temperatures.values[2:(nb.values.variable+1)],all.pdf.P[,H],type="l",
     col=P.col.line,lwd=3,xlim=c(1,max(scale.temperatures.values)),
     main="(c) - P.d.f. of atm. temperature in 2100",cex.main=1.5,cex.axis=1.5,
     cex.lab=1.5,
     xlab="",ylab="",yaxt="no")
lines(scale.temperatures.values[2:(nb.values.variable+1)],all.pdf.Q[,H],
      col=Q.col.line,lwd=3)
abline(v=scale.temperatures.values[which.min(abs(all.cdf.Q[,model_sol$horiz.2100]
                                                 -0.5))],
       col=Q.col.line,lty=3,lwd=2)                                              #median of Q
abline(v=ET.Q[H],col=Q.col.line,lty=1,lwd=1)                                    #mean of Q
abline(v=scale.temperatures.values[which.min(abs(all.cdf.P[,model_sol$horiz.2100]
                                                 -0.5))],
       col=P.col.line,lty=3,lwd=2)                                              #median of P
abline(v=ET.P[H],col=P.col.line,lty=1,lwd=1)                                    #mean of P

legend("topright",
       legend=c("Physical p.d.f.","Risk-adjusted p.d.f.","Mean","Median"),
       lty=c(1,1,1,3),
       col=c(P.col.line,Q.col.line,"black","black"),cex=1.5,
       lwd=c(2,2,1,2),bty = "n")
dev.off()
