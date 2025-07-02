# ==============================================================================
# FIGURE 4. Conditional distribution of future temperatures
# Figure_Tat_P_and_Q_vector_CI.pdf
# ==============================================================================

H <- model_sol$horiz.2100

# Prepare omega vectors (for pricing):
omega_ZCB <- matrix(0,model_sol$n.X)
omega_T.at <- omega_ZCB
omega_T.at[which(model_sol$names.var.X=="T_at")] <- 1
TAT <- which(model_sol$names.var.X=="T_at")

#RCP data (outputs from MAGICC6.0)
RCP_MAGICC <- read.csv("data/RCP_Mat_MAGICC.csv", header=FALSE)
T_RCP30 <- RCP_MAGICC$V6
T_RCP45 <- RCP_MAGICC$V7
T_RCP60 <- RCP_MAGICC$V8
T_RCP85 <- RCP_MAGICC$V9

temp<-read.table("./data/mean_ssp.txt",header=TRUE)
temp_graph<-temp[3:11,]

x <- exp(seq(-5,5,length.out = 1000)) #grid for Fourier

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

# ------------------------------------------------------------------------------
# Plot ----
FILE = paste("/outputs/Figures/Figure_Tat_P_and_Q_vector_CI.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=6, height=6)

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE),
       widths=c(1,1), heights=c(1,1))
par(plt=c(.16,.96,.15,.85))

y.lim <- c(.5,6)
x.lim <- c(2035,2100)

plot(model_sol$vec_date[2:(H+1)],ET.P,
     xlim=x.lim,ylim=y.lim, cex.main=1.5,cex.axis=1.5,cex.lab=1.5,
     col="white",xlab="",ylab="Degrees Celsius",las=1,
     main="(a) - Atmosph. temperature")
grid()

for(i in length(vector.of.CI):1){
  # P
  polygon(c(model_sol$vec_date[2:(H+1)],rev(model_sol$vec_date[2:(H+1)])),
          c(all.CI.P[1,,i],rev(all.CI.P[2,,i])),
          col=P.col,border = NA)
}
lines(model_sol$vec_date[2:(H+1)],
      ET.P,lwd=2,col=P.col.line)


lines(RCP_MAGICC$V1[RCP_MAGICC$V1<2100],
      T_RCP45[RCP_MAGICC$V1<2100],lty=2,lwd=2,col="lightsteelblue4")
lines(RCP_MAGICC$V1[RCP_MAGICC$V1<2100],
      T_RCP60[RCP_MAGICC$V1<2100],lty=4,lwd=2,col="lightsteelblue4")
lines(RCP_MAGICC$V1[RCP_MAGICC$V1<2100],
      T_RCP85[RCP_MAGICC$V1<2100],lty=3,lwd=2,col="lightsteelblue4")
legend("topleft",
       legend=c("RCP4.5","RCP6.0","RCP8.5"),
       lty=c(2,4,3),
       col=c("lightsteelblue4"),
       cex=1.5,
       lwd=2,seg.len = 4,bty = "n")

plot(model_sol$vec_date[2:(H+1)],ET.Q,
     xlim=x.lim,ylim=y.lim,las=1,
     col="white",xlab="",ylab="Degrees Celsius", cex.main=1.5,cex.axis=1.5,
     cex.lab=1.5,
     main="(b) - Atmosph. temperature (risk-adjusted)")
grid()

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


par(plt=c(.08,.98,.15,.85))

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
