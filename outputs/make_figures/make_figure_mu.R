# ==============================================================================
# Figure showing mitigation path (mu_t)
# ==============================================================================
# Check sensitivity to assumptions regarding mu_t:
# (i) parametric function + (ii) time consistency
# ==============================================================================


x.lim <- c(model_sol$vec_date[2],2200)

dice_data           <- read.csv("data/mu.csv",sep=";",header = F)
row.names(dice_data)<- dice_data[,1]
dice_data           <- dice_data[-1]

mu_dice                    <-matrix(0,dim(dice_data)[2],1)
mu_dice[1:length(mu_dice)] <-apply(dice_data[2,],2,function(x)min(x,1))

#Plot
FILE = paste("/outputs/Figures/Figure_Mitigation_comparison.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=9,width=7, height=5)

par(mfrow=c(2,2))
par(plt=c(.15,.95,.15,.8))

plot(model_sol$vec_date[2:length(model_sol$vec_date)],
     model_sol$mu[1:(length(model_sol$vec_date)-1)],
     type="l", xlab="Year",ylab="",lwd=2,
     ylim=c(0,1),las=1,
     xlim=x.lim,las=1,main="(a) Mitigation rate, comparison with DICE")
lines(model_sol$vec_date[2:length(model_sol$vec_date)],
      mu_dice[2:length(model_sol$vec_date)],lwd=2,lty=2,
      col="grey")


# (i) Check sensitivity to parametric function ---------------------------------

t <- 0
Xt <- model_sol$X

RES.2param <- res.optim(model_sol,
                        model_sol$theta0,
                        Tend = model_sol$Tmax - t,
                        X = Xt)

# use of a logit function to code mu_t:
theta <- .5 + .99*(model_sol$mu[(t+1):model_sol$Tmax] - .5)
theta <- log(theta/(1-theta))
RES.Tmax.param <- res.optim(model_sol,
                            theta,
                            Tend = model_sol$Tmax - t,
                            X = Xt)

legend("bottomright",
       legend=c("Parametric function","Non-parametric function","DICE"),
       lty=c(1,NaN,2),
       col=c("black","black","grey"),
       lwd=c(2,1,2),
       pch=c(NaN,3,NaN),cex=1)

points(model_sol$vec_date[2:length(model_sol$vec_date)],
       mu.function(model_sol,
                  theta = RES.Tmax.param$par,
                  t.ini = t)[1:(length(model_sol$vec_date)-1)],
       col="black",pch=3)


# (ii) Check sensitivity to no-update assumption -------------------------------

EV <- EV.fct(model_sol)

max.t <- 40

# Run parallel computations:
cl <- makeCluster(number.of.cores)
registerDoParallel(cl)
save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))
clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))

# --- State vector = avg
all_mu <- foreach(t = 1:max.t, 
                   .combine=rbind) %dopar% {
                     
                     Xt <- EV$EXh[[t]]
                     
                     Xt <- Xt
                     
                     # Re-optimize future trajectory of mu_t's:
                     RES.2param <- res.optim(model_sol,
                                             model_sol$theta0,
                                             Tend = model_sol$Tmax - t,
                                             X = Xt)
                     
                     # Compute resulting mu_t:
                     mu <- mu.function(model_sol,
                                       theta = RES.2param$par,
                                       t.ini = t)
                     mu
                   }
plot(model_sol$vec_date[2:length(model_sol$vec_date)],
     mu.function(model_sol,theta = RES.2param$par,
                 t.ini = t)[1:(length(model_sol$vec_date)-1)],
     type="l",col="white",lwd=3,
     ylim=c(0,1),las=1,xlab="years",ylab="",
     main="(b) Re-optimize - State vector = avg")
for(t in 1:max.t){
  lines(model_sol$vec_date[2:length(model_sol$vec_date)],
        all_mu[t,1:(length(model_sol$vec_date)-1)],col="black")
}
lines(model_sol$vec_date[2:length(model_sol$vec_date)],
      mu.function(model_sol,
                 theta = RES.2param$par,
                 t.ini = t)[1:(length(model_sol$vec_date)-1)],
      type="l",col="dark grey",lty=3,lwd=4)
legend("bottomright",
       legend=c("Baseline case","Re-optimized rates"),
       lty=c(3,1),
       col=c("grey","black"),
       lwd=c(4,1),cex=1)

# --- State vector = avg + 1std dev
all_mu <- foreach(t = 1:max.t, 
                  .combine=rbind) %dopar% {
                    
                    Xt <- EV$EXh[[t]]
                    
                    Xt <- Xt + sqrt(diag(EV$CovX[[t]]))
                    
                    # Re-optimize future trajectory of mu_t's:
                    RES.2param <- res.optim(model_sol,
                                            model_sol$theta0,
                                            Tend = model_sol$Tmax - t,
                                            X = Xt)
                    
                    # Compute resulting mu_t:
                    mu <- mu.function(model_sol,
                                      theta = RES.2param$par,
                                      t.ini = t)
                    mu
                  }
plot(model_sol$vec_date[2:length(model_sol$vec_date)],
     mu.function(model_sol,theta = RES.2param$par,
                 t.ini = t)[1:(length(model_sol$vec_date)-1)],
     type="l",col="white",lwd=3,
     ylim=c(0,1),las=1,xlab="years",ylab="",
     main="(c) Re-optimize - State vector = avg + 1 std.dev.")
for(t in 1:max.t){
  lines(model_sol$vec_date[2:length(model_sol$vec_date)],
        all_mu[t,1:(length(model_sol$vec_date)-1)],col="black")
}
lines(model_sol$vec_date[2:length(model_sol$vec_date)],
      mu.function(model_sol,
                  theta = RES.2param$par,
                  t.ini = t)[1:(length(model_sol$vec_date)-1)],
      type="l",col="dark grey",lty=3,lwd=4)
legend("bottomright",
       legend=c("Baseline case","Re-optimized rates"),
       lty=c(3,1),
       col=c("grey","black"),
       lwd=c(4,1),cex=1)

# --- State vector = avg - 1std dev
all_mu <- foreach(t = 1:max.t, 
                  .combine=rbind) %dopar% {
                    
                    Xt <- EV$EXh[[t]]
                    
                    Xt <- Xt - sqrt(diag(EV$CovX[[t]]))
                    
                    # Re-optimize future trajectory of mu_t's:
                    RES.2param <- res.optim(model_sol,
                                            model_sol$theta0,
                                            Tend = model_sol$Tmax - t,
                                            X = Xt)
                    
                    # Compute resulting mu_t:
                    mu <- mu.function(model_sol,
                                      theta = RES.2param$par,
                                      t.ini = t)
                    mu
                  }
plot(model_sol$vec_date[2:length(model_sol$vec_date)],
     mu.function(model_sol,theta = RES.2param$par,
                 t.ini = t)[1:(length(model_sol$vec_date)-1)],
     type="l",col="white",lwd=3,
     ylim=c(0,1),las=1,xlab="years",ylab="",
     main="(d) Re-optimize - State vector = avg - 1 std.dev.")
for(t in 1:max.t){
  lines(model_sol$vec_date[2:length(model_sol$vec_date)],
        all_mu[t,1:(length(model_sol$vec_date)-1)],col="black")
}
lines(model_sol$vec_date[2:length(model_sol$vec_date)],
      mu.function(model_sol,
                  theta = RES.2param$par,
                  t.ini = t)[1:(length(model_sol$vec_date)-1)],
      type="l",col="dark grey",lty=3,lwd=4)
legend("bottomright",
       legend=c("Baseline case","Re-optimized rates"),
       lty=c(3,1),
       col=c("grey","black"),
       lwd=c(4,1),cex=1)

stopCluster(cl)
file.remove("outputs/toto.Rdata")

dev.off()
