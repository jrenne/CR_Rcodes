# ==============================================================================
# FIGURE S.2. Mitigation rate
# Figure_Mitigation_comparison.pdf
# ==============================================================================


# Check sensitivity to assumptions regarding mu_t:
# (i) parametric function + (ii) time consistency

x.lim <- c(model_sol$vec_date[2],2150)

mu_DICE <- read.csv("data/mu_DICE.csv")

indic_Mat     <- which(model_sol$names.var.X=="M_at")
indic_Mup     <- which(model_sol$names.var.X=="M_up")
indic_Mlo     <- which(model_sol$names.var.X=="M_lo")
indic_N       <- which(model_sol$names.var.X=="N")
indic_forc    <- which(model_sol$names.var.X=="Forc")
indic_T_at    <- which(model_sol$names.var.X=="T_at")
indic_y_tilde <- which(model_sol$names.var.X=="y_tilde")
indics <- c(indic_Mat,
            indic_Mup,
            indic_Mlo,
            indic_T_at,
            indic_y_tilde) # changes in Xt will be applied to these entries
indics <- 1:model_sol$n.X

# to compute temperature risk premium:
omega_ZCB <- matrix(0,model_sol$n.X,1)
omega_T.at <- omega_ZCB
omega_T.at[indic_T_at] <- 1

model_sol$theta0 <- model_sol$theta.opt

# ------------------------------------------------------------------------------
# Plot----
FILE = paste("/outputs/Figures/Figure_Mitigation_comparison.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=9,width=7, height=5)

par(mfrow=c(2,2))
par(plt=c(.15,.95,.15,.8))

plot(model_sol$vec_date[2:length(model_sol$vec_date)],
     model_sol$mu[2:length(model_sol$vec_date)],
     type="l", xlab="Year",ylab="",lwd=2,
     ylim=c(0,1.1),las=1,
     xlim=x.lim,
     las=1,main="(a) Mitigation rate, comparison with DICE")

lines(mu_DICE$X.Year,
      mu_DICE$Emissions.Control.Rate,lwd=2,lty=2,
      col="grey")

# (i) Check sensitivity to parametric function ---------------------------------

t <- 0
Xt <- model_sol$X

# Baseline approach:
RES.2param <- res.optim(model_sol,
                        model_sol$theta0,
                        Tend = model_sol$Tmax - t,
                        X = Xt)

# Non-parametric function:
theta <- .5 + .99*(model_sol$mu[(t+1):model_sol$Tmax] - .5)
theta <- log(theta/(1-theta)) # starting values
RES.Tmax.param <- res.optim(model_sol,
                            theta,
                            Tend = model_sol$Tmax - t,
                            X = Xt)

legend("bottomright",
       legend=c("Parametric function","Non-parametric function","DICE (2023)"),
       lty=c(1,NaN,2),
       col=c("black","black","grey"),
       lwd=c(2,1,2),
       pch=c(NaN,3,NaN),cex=1)

points(model_sol$vec_date[2:length(model_sol$vec_date)],
       mu.function(model_sol,
                   theta = RES.Tmax.param$par,
                   t.ini = t)[2:length(model_sol$vec_date)],
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
                    Xt[indics] <- Xt[indics] + 0*sqrt(diag(EV$CovX[[t]]))[indics]

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
                 t.ini = t)[2:length(model_sol$vec_date)],
     type="l",col="white",lwd=3,
     ylim=c(0,1.1),
     xlim=x.lim,
     las=1,xlab="years",ylab="",
     main="(b) Re-optimize - State vector = avg")
for(t in 1:max.t){
  lines(model_sol$vec_date[2:length(model_sol$vec_date)],
        all_mu[t,2:length(model_sol$vec_date)],col="black")
}
lines(model_sol$vec_date[2:length(model_sol$vec_date)],
      mu.function(model_sol,
                  theta = RES.2param$par,
                  t.ini = t)[2:length(model_sol$vec_date)],
      type="l",col="dark grey",lty=3,lwd=4)
legend("bottomright",
       legend=c("Baseline case","Re-optimized rates"),
       lty=c(3,1),
       col=c("grey","black"),
       lwd=c(4,1),cex=1)

# --- State vector = avg + 2std dev
all_mu <- foreach(t = 1:max.t, 
                  .combine=rbind) %dopar% {
                    
                    Xt <- EV$EXh[[t]]
                    Xt[indics] <- Xt[indics] - 2*sqrt(diag(EV$CovX[[t]]))[indics]

                    # Re-optimize future trajectory of mu_t's:
                    RES.2param <- res.optim(model_sol,
                                            model_sol$theta0,
                                            Tend = model_sol$Tmax - t,
                                            X = Xt)
                    
                    # Compute resulting mu_t:
                    mu <- mu.function(model_sol,
                                      theta = RES.2param$par,
                                      t.ini = t)
                    
                    # Recover new utility specification (for SCC analysis):-----
                    u_specif <- 
                      utility.optim(model_sol,
                                    theta = RES.2param$par,
                                    Tend = model_sol$Tmax - t,
                                    X = Xt,
                                    indic_returns_mu_u.1 = TRUE)
                    u_specif_baseline <- 
                      utility.optim(model_sol,
                                    theta = model_sol$theta.opt,
                                    Tend = model_sol$Tmax - t,
                                    X = Xt,
                                    indic_returns_mu_u.1 = TRUE)
                    # Temperature RP: ------
                    # with new model (after updating of mu):
                    HHH <- 10 # maturity in model periods
                    model_sol_new <- model_sol
                    model_sol_new <- model_solve(model_sol_new,indic_mitig = FALSE,
                                                 mu.chosen = mu.function(model_sol,
                                                                         theta = RES.2param$par))
                    model_sol_new_P <- model_sol_new
                    # Kill pricing of risks in model_sol_P:
                    model_sol_new_P$pi <- lapply(model_sol_new_P$pi,function(x){0*x})
                    model_sol_new_P$eta0 <- 0*model_sol_new_P$eta0
                    model_sol_new_P$eta1 <- lapply(model_sol_new_P$eta1,function(x){0*x})
                    ET.P_new  <- varphi.tilde(model_sol_new_P,omega_T.at,HHH,
                                              X=Xt,t=t)[[1]]
                    ZCB <- varphi(model_sol_new,omega_ZCB,HHH,X=Xt,t=t)
                    ET.Q_new  <- varphi.tilde(model_sol_new,omega_T.at,HHH,
                                              X=Xt,t=t)[[1]]/ZCB$P.t
                  
                    # compute TRP for date t with baseline model:
                    # with baseline model (without updating of mu):
                    model_sol_P <- model_sol
                    # Kill pricing of risks in model_sol_P:
                    model_sol_P$pi <- lapply(model_sol_P$pi,function(x){0*x})
                    model_sol_P$eta0 <- 0*model_sol_P$eta0
                    model_sol_P$eta1 <- lapply(model_sol_P$eta1,function(x){0*x})
                    ET.P  <- varphi.tilde(model_sol_P,omega_T.at,HHH,
                                          X=Xt,t=t)[[1]]
                    ZCB <- varphi(model_sol,omega_ZCB,HHH,X=Xt,t=t)
                    ET.Q  <- varphi.tilde(model_sol,omega_T.at,HHH,
                                          X=Xt,t=t)[[1]]/ZCB$P.t
                    # plot(ET.P,type="l")
                    # lines(ET.Q,lty=2)
                    # lines(ET.P_new,col="red")
                    # lines(ET.Q_new,col="red",lty=2)
                    # ----------------------------------------------------------
                    
                    # The last two entries are to check the effect of updating on SCC:
                    c(mu,
                      u_specif[[2]][indic_Mat],
                      u_specif_baseline[[2]][indic_Mat],
                      ET.Q[HHH] - ET.P[HHH],ET.Q_new[HHH] - ET.P_new[HHH],
                      ET.Q[HHH],ET.P[HHH],ET.Q_new[HHH],ET.P_new[HHH])
                  }
plot(model_sol$vec_date[2:length(model_sol$vec_date)],
     mu.function(model_sol,theta = RES.2param$par,
                 t.ini = t)[2:length(model_sol$vec_date)],
     type="l",col="white",lwd=3,
     ylim=c(0,1.1),
     xlim=x.lim,
     las=1,xlab="years",ylab="",
     main="(c) Re-optimize - State vector = avg + 2 std.dev.")
for(t in 1:max.t){
  lines(model_sol$vec_date[2:length(model_sol$vec_date)],
        all_mu[t,2:length(model_sol$vec_date)],col="black")
}
lines(model_sol$vec_date[2:length(model_sol$vec_date)],
      mu.function(model_sol,
                  theta = RES.2param$par,
                  t.ini = t)[2:length(model_sol$vec_date)],
      type="l",col="dark grey",lty=3,lwd=4)
legend("bottomright",
       legend=c("Baseline case","Re-optimized rates"),
       lty=c(3,1),
       col=c("grey","black"),
       lwd=c(4,1),cex=1)

# compute change in SCC wrt basline:
chge_SCC <- (all_mu[,model_sol$Tmax+2] - all_mu[,model_sol$Tmax+1])/
  all_mu[,model_sol$Tmax+1]
chge_TRP <- (all_mu[,model_sol$Tmax+4] - all_mu[,model_sol$Tmax+3])/
  all_mu[,model_sol$Tmax+3]
# print(paste("Average absolute percent change in SCC: ",round(mean(abs(chge_SCC)),5),sep=""))
# print(paste("Average absolute percent change in TRP: ",round(mean(abs(chge_TRP)),5),sep=""))

# --- State vector = avg - 2std dev
all_mu <- foreach(t = 1:max.t, 
                  .combine=rbind) %dopar% {
                    
                    Xt <- EV$EXh[[t]]
                    Xt[indics] <- Xt[indics] - 2*sqrt(diag(EV$CovX[[t]]))[indics]

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
                 t.ini = t)[2:length(model_sol$vec_date)],
     type="l",col="white",lwd=3,
     ylim=c(0,1.1),
     xlim=x.lim,
     las=1,xlab="years",ylab="",
     main="(d) Re-optimize - State vector = avg - 2 std.dev.")
for(t in 1:max.t){
  lines(model_sol$vec_date[2:length(model_sol$vec_date)],
        all_mu[t,2:length(model_sol$vec_date)],col="black")
}
lines(model_sol$vec_date[2:length(model_sol$vec_date)],
      mu.function(model_sol,
                  theta = RES.2param$par,
                  t.ini = t)[2:length(model_sol$vec_date)],
      type="l",col="dark grey",lty=3,lwd=4)
legend("bottomright",
       legend=c("Baseline case","Re-optimized rates"),
       lty=c(3,1),
       col=c("grey","black"),
       lwd=c(4,1),cex=1)

dev.off()

stopCluster(cl)
file.remove("outputs/toto.Rdata")
