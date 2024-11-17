# ==============================================================================
# Term Structure of Yields, BEIR
# ==============================================================================

H <- 100

# Prepare omega vectors (for pricing):
omega_ZCB <- matrix(0,model_sol$n.X)
omega_T.at <- omega_ZCB
omega_T.at[which(model_sol$names.var.X=="T_at")] <- 1

mu.pi.0     <- 0.02*model$tstep
all.mu.pi.D <- c(0.01,1,5)

colors.muD <- sequential_hcl(length(all.mu.pi.D)+2, "YlOrRd")
colors.muD <- colors.muD[length(all.mu.pi.D):1]

indic_D  <- which(model$names.var.X=="D")
indic_Pi <- which(model$names.var.X=="Cum_dc")

# Compute term structure of real yields:
prices.ZCRF.bonds <- varphi(model_sol,omega.varphi = omega_ZCB,H=H)

# Determine expectations of damages (to adjust ex-damages inflation):
EV <- EV.fct(model_sol)
ED <- EV$EX$D

all.nom.rates <- NULL
all.infl.exp  <- NULL
all.BEIR      <- NULL


names4legend <- c("Real yields")
for(k in 1:length(all.mu.pi.D)){
  eval(parse(text = gsub(" "," ",
                         paste("names4legend <- c(names4legend,expression(paste('Nominal yields (',mu[pi*','*D],' = ',",all.mu.pi.D[k],",')',sep='')))")
  )))
}

panel.titles <- NULL
for(k in 1:length(all.mu.pi.D)){
  eval(parse(text = gsub(" "," ",
                         paste("panel.titles <- c(panel.titles,expression(paste('Panel (",letters[k+1],") BEIR with ',mu[pi*','*D],' = ',",all.mu.pi.D[k],",sep='')))",sep="")
  )))
}


# Compute term structures of nominal yields
for(k in 1:length(all.mu.pi.D)){
  mu_pi_0_t <- mu.pi.0 - all.mu.pi.D[k] * ED
  if(H>length(mu_pi_0_t)){
    mu_pi_0_t <- c(mu_pi_0_t,rep(mu_pi_0_t[length(mu_pi_0_t)],H-length(mu_pi_0_t)))
  }
  mu_PI <- list(muprice_0 = mu_pi_0_t,
                muprice_1 = matrix(0,model_sol$n.X,1))
  mu_PI$muprice_1[indic_D] <- all.mu.pi.D[k]
  model_sol <- update.model_sol.4.mu_altern(model_sol,
                                            mu_altern=mu_PI)
    
  omega_PI <- matrix(0,model_sol$n.X,1)
  omega_PI[indic_Pi] <- 1
  
  prices.nom.bonds <- varphi(model_sol,omega.varphi = - omega_PI,H=H)
  nom.rates        <- prices.nom.bonds$r.t - 
    100/((1:H)*model$tstep)*c(t(mu_PI$muprice_1)%*%model_sol$X)
  
  #Model-implied EV
  EV<-EV.fct(model_sol,H)
  # Get matrix of forecasts:
  EXh <- matrix(unlist(EV$EXh),model_sol$n.X,H)
  expected.annualized.inflation <- c(100*(mu_PI$muprice_0 +
                                          matrix(mu_PI$muprice_1,nrow = 1) %*% EXh)/model$tstep)
  
  all.nom.rates <- cbind(all.nom.rates,nom.rates)
  all.infl.exp  <- cbind(all.infl.exp,expected.annualized.inflation)
  all.BEIR      <- cbind(all.BEIR,nom.rates - prices.ZCRF.bonds$r.t)
}


# Plot ----
FILE = "/outputs/Figures/Figure_BreakEvenInflation.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=11, width=7, height=6)
par(plt=c(.1,.95,.25,.85))


par(mfrow=c(2,1))
plot(model$tstep*(1:H),prices.ZCRF.bonds$r.t,
     ylim=c(1.5,trunc(max(all.nom.rates,na.rm = TRUE))+1),
     type="l",lwd=3,
     xlab="maturity (in years)",ylab="yield-to-maturity (in percent)",
     cex.lab=.9,cex.main=.9,cex.axis=.9,
     lty=3,las=1,
     main=expression(paste("Panel (a) Term structures of nominal and real interest rates",sep="")))
for(k in 1:length(all.mu.pi.D)){
  lines(model$tstep*(1:H),all.nom.rates[,k],
        lwd=2,col=colors.muD[k])
}
grid()

par(mfrow=c(2,length(all.mu.pi.D)))

for(k in 1:length(all.mu.pi.D)){
  
  par(plt=c(.25,.95,.25,.85))
  
  par(mfg=c(2,k))
  
  main.title <- expression()
    
  plot(model$tstep*(1:H),all.infl.exp[,k],
       ylim=c(min(all.infl.exp,na.rm = TRUE),
              max(all.BEIR,na.rm = TRUE)),
       type="l",lwd=2,lty=2,las=1,
       cex.lab=1.2,cex.main=1.2,cex.axis=1.2,
       xlab="maturity (in years)",ylab="inflation (in percent)",
       main = panel.titles[k])
  lines(model$tstep*(1:H),all.BEIR[,k],
        lwd=2,col=colors.muD[k])
}

dev.off()





