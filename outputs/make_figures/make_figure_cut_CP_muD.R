# Sensitivity of Temperature risk premium to mu_D:

H <- model_sol$horiz.2100

# Prepare omega vectors (for pricing):
omega_ZCB <- matrix(0,model_sol$n.X)
omega_T.at <- omega_ZCB
omega_T.at[which(model_sol$names.var.X=="T_at")] <- 1

values.of.mu_D <- seq(0.6*model_sol$parameters$mu_D,
                      2*model_sol$parameters$mu_D,length.out = 24)

cl <- makeCluster(number.of.cores)
registerDoParallel(cl)

save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))

clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))

all.Ts <- foreach(i = 1:length(values.of.mu_D), .combine=rbind) %dopar% {
  
  model_i <- model_sol
  
  model_i$parameters$mu_D <- values.of.mu_D[i]
  model_i_sol <- model_solve(model_i,theta0)
  
  EV.i       <- EV.fct(model_i_sol,h=H)
  T.P        <- EV.i$EX$T_at[H]
  varphi.i   <- varphi.tilde(model_i_sol,omega_T.at,H)[[1]]/
    varphi(model_i_sol,omega_ZCB,H)[[3]]
  T.Q        <- varphi.i[H]
  c(T.P,T.Q)
}

stopCluster(cl)
file.remove("outputs/toto.Rdata")


T.P       <- matrix(all.Ts[,1],length(values.of.mu_D),1)
T.Q       <- matrix(all.Ts[,2],length(values.of.mu_D),1)

#Plots
#1
FILE = paste("/outputs/Figures/Figure_cut_CP.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=7, height=3.5)

par(mfrow=c(1,1))
par(plt=c(.1,.95,.15,.95))

# Refine grid (interpolation):
values.of.mu_D.fine <- seq(values.of.mu_D[1],
                           tail(values.of.mu_D,1),length.out = 50)

# # Use splines:
# spline.T.P<-splinefun(values.of.mu_D,T.P,method="natural")
# spline.T.Q<-splinefun(values.of.mu_D,T.Q,method="natural")
# T.P.fit <- spline.T.P(values.of.mu_D.fine)
# T.Q.fit <- spline.T.Q(values.of.mu_D.fine)

# Other type of splines:
# spline.T.P <- smooth.spline(values.of.mu_D,T.P,df=4)
# T.P.fit    <- predict(spline.T.P,values.of.mu_D.fine)$y
# spline.T.Q <- smooth.spline(values.of.mu_D,T.Q,df=4)
# T.Q.fit    <- predict(spline.T.Q,values.of.mu_D.fine)$y

values.of.mu_D.fine <- values.of.mu_D
T.P.fit <- T.P
T.Q.fit <- T.Q

plot(values.of.mu_D.fine,T.P.fit,type="l",
     col=P.col.line,lwd=2,las=1,
     cex.main=1.5,cex.axis=1.5,cex.lab=1.5,
     main="",
     xlab=expression(paste("Disaster uncertainty ",mu[D],sep="")),
     ylab=expression(paste("Temperature Anomaly ",T[at],sep="")),
     ylim=c(min(T.P.fit,T.Q.fit),
            min(5,max(T.P.fit,T.Q.fit))),lty=1)
# lines(values.of.mu_D.fine,T.Q.fit,
#       col=Q.col.line,lwd=2,lty=1)
points(values.of.mu_D,T.Q,
      col=Q.col.line,lwd=2,pch=16,cex=1.3)
abline(v=model_sol$parameters$mu_D,col="grey",lty=3,lwd=2)
legend("topleft",
       legend=c(expression(paste("Expected ",T[at]," in 2100",sep="")),
                expression(paste("Swaps Price ",T^S," in 2100",sep=""))
       ),
       lty=c(1,NaN),
       pch=c(NaN,16),
       lwd=c(2,2),
       col=c(P.col.line,Q.col.line),
       bty = "n",cex=1.4,seg.len = 3)
dev.off()

