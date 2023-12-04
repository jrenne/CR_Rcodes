# ==============================================================================
# Confidence interval of temperature risk premiums
# ==============================================================================

max.chge <- .5 # maximum absolute change in target

Nb.draws <- 100

H <- model_sol$horiz.2100

gammas <- c(model$parameters$gamma,2)
deltas <- c(model$parameters$delta,(1-.01)^model$tstep)

rhos <- 1 - deltas^(1/model$tstep)

all.gamma <- gammas %x% c(1,1)
all.delta <- c(1,1) %x% deltas

# Prepare omega vectors (for pricing):
omega_ZCB  <- matrix(0,model_sol$n.X)
omega_T.at <- matrix(0,model_sol$n.X)
omega_T.at[which(model_sol$names.var.X=="T_at")] <- 1
omega_H <- matrix(0,model_sol$n.X)
omega_H[which(model_sol$names.var.X=="H")] <- 1

all.SCC <- matrix(NaN,Nb.draws,length(all.gamma))
all.TRP <- matrix(NaN,Nb.draws,length(all.gamma))
all.HRP <- matrix(NaN,Nb.draws,length(all.gamma))
all.LTR <- matrix(NaN,Nb.draws,length(all.gamma))


for(iii in 1:length(all.gamma)){
  
  print(paste("Progression: working on case ",iii," (out of ",
              length(all.gamma),")",sep=""))
  
  cl <- makeCluster(number.of.cores)
  registerDoParallel(cl)
  
  save.image("outputs/toto.Rdata")
  clusterEvalQ(cl,load("outputs/toto.Rdata"))
  
  clusterEvalQ(cl,library(MASS))
  clusterEvalQ(cl,library(expm))
  
  all.res <- foreach(i = 1:Nb.draws, .combine=rbind) %dopar% {
    
    set.seed(i)
    
    model_sol_new <- model_sol
    model_sol_new$parameters$delta <- all.delta[iii]
    model_sol_new$parameters$gamma <- all.gamma[iii]
    
    # Damage targets:
    model_sol_new$target_vector["ECumD2"] <- 1 - (1-model_sol$target_vector["ECumD2"]) *
      (1 + max.chge*runif(1,min=-1,max=1))
    model_sol_new$target_vector["ECumD4"] <- model_sol_new$target_vector["ECumD2"] -
      (model_sol$target_vector["ECumD2"] - model_sol$target_vector["ECumD4"]) *
      (1 + max.chge*runif(1,min=-1,max=1))
    model_sol_new$target_vector["stdCumD4"] <- model_sol$target_vector["stdCumD4"] *
      (1 + max.chge*runif(1,min=-1,max=1))
    
    # Permafrost-related targets:
    model_sol_new$target_vector["ECumN2"] <- model_sol$target_vector["ECumN2"] *
      (1 + max.chge*runif(1,min=-1,max=1))
    model_sol_new$target_vector["ECumN4"] <- model_sol_new$target_vector["ECumN2"] +
      (model_sol$target_vector["ECumN4"] - model_sol$target_vector["ECumN2"]) *
      (1 + max.chge*runif(1,min=-1,max=1))
    model_sol_new$target_vector["stdCumN4"] <- model_sol$target_vector["stdCumN4"] *
      (1 + max.chge*runif(1,min=-1,max=1))
    model_sol_new$target_vector["ECumNinf"] <- model_sol$target_vector["ECumNinf"] *
      (1 + max.chge*runif(1,min=-1,max=1))
    
    # SLR targets:
    model_sol_new$target_vector["EH2"] <- model_sol$target_vector["EH2"] *
      (1 + max.chge*runif(1,min=-1,max=1))
    model_sol_new$target_vector["EH4"] <- model_sol_new$target_vector["EH2"] +
      (model_sol$target_vector["EH4"] - model_sol$target_vector["EH2"]) *
      (1 + max.chge*runif(1,min=-1,max=1))
    model_sol_new$target_vector["stdH4"] <- model_sol$target_vector["stdH4"] *
      (1 + max.chge*runif(1,min=-1,max=1))
    
    # Consumption moment targets:
    model_sol_new$target_vector["mu_c0"] <- model_sol$target_vector["mu_c0"] *
      (1 + max.chge*runif(1,min=-1,max=1))
    model_sol_new$target_vector["sigma_c0"] <- model_sol$target_vector["sigma_c0"] *
      (1 + max.chge*runif(1,min=-1,max=1))
    
    # Temperature dynamics uncertainty:
    mu_T <- model_sol_new$parameters$mu_T * (1 + max.chge*runif(1,min=-1,max=1))
    model_sol_new$parameters$mu_T <- mu_T
    
    targets_and_mu <- c(model_sol_new$target_vector,mu_T)
    save(targets_and_mu,
         file=paste("outputs/results/targets_and_mu",i,".Rdat",sep=""))
    
    model_sol_new <- solveParam4D(model_sol_new)
    model_sol_new <- solveParam4H(model_sol_new)
    model_sol_new <- solveParam4N(model_sol_new)
    model_sol_new <- solveParam4c(model_sol_new)
    
    model_sol_new <- model_solve(model_sol_new,indic_CRRA = FALSE)
    
    SCC <- scc.fct(model_sol_new,h=0)
    
    EV    <- EV.fct(model_sol_new,h=H)
    
    ZCB <- varphi(model_sol_new,omega_ZCB,H)
    
    ET.P  <- EV$EX$T_at[1:H]
    ET.Q  <- 
      varphi.tilde(model_sol_new,omega_T.at,H)[[1]]/ZCB$P.t
    
    EH.P  <- EV$EX$H[1:H]
    EH.Q  <- 
      varphi.tilde(model_sol_new,omega_H,H)[[1]]/ZCB$P.t
    
    T.RP <- c(ET.Q - ET.P)
    H.RP <- c(EH.Q - EH.P)
    
    c(SCC,T.RP,H.RP,ZCB$r.t)
  }
  
  stopCluster(cl)
  file.remove("outputs/toto.Rdata")
  
  all.SCC[,iii] <- all.res[,1]
  all.TRP[,iii] <- all.res[,1 + 0*H + H]
  all.HRP[,iii] <- all.res[,1 + 1*H + H]
  all.LTR[,iii] <- all.res[,1 + 2*H + H]
  
}



label.gamma <- NULL
label.delta <- NULL
for(j in 1:length(gammas)){
  eval(parse(text = gsub(" "," ",
                         paste("label <- expression(paste(gamma,' = ',",
                               gammas[j],",sep=''))",sep="")
  )))
  label.gamma <- c(label.gamma,label)
  # eval(parse(text = gsub(" "," ",
  #                        paste("label <- expression(paste(delta,' = (1 - ",
  #                              rhos[j],")^(5)',sep=''))",sep="")
  # )))
  eval(parse(text = gsub(" "," ",
                         paste("label <- expression(paste('Discount factor ',delta,' = ',(1-",
                               100*rhos[j],"/100)^5,sep=''))",sep="")
  )))
  label.delta <- c(label.delta,label)
}


FILE = paste("/outputs/Figures/Figure_IntConf_RP.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=11,width=9, height=6)

par(mfrow=c(2,2))
plot.new()
title("(a) Scocial Cost of Carbon (in U.S. $)")
plot.new()
title("(b) Temperature risk premium (in °C)")
plot.new()
title("(c) SLR risk premium (in meters)")
plot.new()
title("(d) Long-term rate (in percent)")

par(mfrow=c(2,4))
par(mgp=c(6,2,0)) 
par(plt=c(.2,.95,.2,.7))

# SCC
par(mfg=c(1,1))
boxplot(all.SCC[,c(1,3)],outline = FALSE,
        names=c(label.gamma),
        main=label.delta[1],las=1)
abline(h=0,col="grey",lty=3)
par(mfg=c(1,2))
boxplot(all.SCC[,c(2,4)],outline = FALSE,
        names=c(label.gamma),
        main=label.delta[2],las=1)
abline(h=0,col="grey",lty=3)
# TRP
par(mfg=c(1,3))
boxplot(all.TRP[,c(1,3)],outline = FALSE,
        names=c(label.gamma),
        main=label.delta[1],las=1)
abline(h=0,col="grey",lty=3)
par(mfg=c(1,4))
boxplot(all.TRP[,c(2,4)],outline = FALSE,
        names=c(label.gamma),
        main=label.delta[2],las=1)
abline(h=0,col="grey",lty=3)
# HRP
par(mfg=c(2,1))
boxplot(all.HRP[,c(1,3)],outline = FALSE,
        names=c(label.gamma),
        main=label.delta[1],las=1)
abline(h=0,col="grey",lty=3)
par(mfg=c(2,2))
boxplot(all.HRP[,c(2,4)],outline = FALSE,
        names=c(label.gamma),
        main=label.delta[2],las=1)
abline(h=0,col="grey",lty=3)
# LTR
par(mfg=c(2,3))
boxplot(all.LTR[,c(1,3)],outline = FALSE,
        names=c(label.gamma),
        main=label.delta[1],las=1)
abline(h=0,col="grey",lty=3)
par(mfg=c(2,4))
boxplot(all.LTR[,c(2,4)],outline = FALSE,
        names=c(label.gamma),
        main=label.delta[2],las=1)
abline(h=0,col="grey",lty=3)

dev.off()

# all.SCC.noNa <- all.SCC[!is.na(all.SCC)]
# mean.SCC <- mean(all.SCC.noNa)
# var.SCC <- var(all.SCC.noNa)
# xs.kurtosis.SCC <- mean((all.SCC.noNa - mean.SCC)^4)/sd(all.SCC.noNa)^4 - 3
# 
# sigma <- seq(0,1,by=.001)
# xs.kurtosis <- exp(4*sigma^2) + 2*exp(3*sigma^2) + 3*exp(2*sigma^2) - 6
# indic.sol <- which((xs.kurtosis - xs.kurtosis.SCC)^2 == min((xs.kurtosis - xs.kurtosis.SCC)^2))
# sigma <- sigma[indic.sol]
# mu    <- (log(var.SCC/(exp(sigma^2)-1)) - sigma^2)/2
# mean.lognorm <- exp(mu + sigma^2/2)
# 
# x <- seq(0,1000,by=1)
# pdf <- dlnorm(x,meanlog=mu,sdlog=sigma)
# plot(x + mean.SCC - mean.lognorm,pdf,type="l")
# abline(v=mean.SCC,col="red")
# lines(density(all.SCC.noNa),lty=2)