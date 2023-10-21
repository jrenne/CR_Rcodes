
ylim.returns <- 100*c(.03,.05)
ylim.PD      <- c(log(.0001),log(.25))
ylim.ratioPD <- c(0,35)
ylim.ratioRates <- c(2.5,4)
xlim         <- c(.2,1.2)
#xlim         <- c(.25,1.5)

# Values of debt to repay:
A_bar <- seq(.01,max(1,max(xlim)),length.out=50)

horizon <- 4 # in number of periods

# For Fourier-transform computation:
x <- exp(seq(-5,5,length.out = 10000)) # grid for Proposition 8 (Fourier)

# Returns of risk-free strategy:
omega.ZC <- matrix(0,model_sol$n.X,1)
prices.ZCRF.bonds   <- varphi(model_sol,
                              omega.varphi = omega.ZC,
                              H = horizon)
RF.strategy <- 1/prices.ZCRF.bonds$P.t


# panel.titles <- NULL
# for(k in 1:length(vector.of.muAD)){
#   eval(parse(text = gsub(" "," ",
#                          paste("panel.titles <- c(panel.titles,expression(paste('Panel (",
#                                letters[k],") ',mu[A*','*H],' = ',",
#                                vector.of.muAD[k],",sep='')))",sep="")
#   )))
# }
legend.muAD <- NULL
for(k in 1:length(vector.of.muAD)){
  eval(parse(text = gsub(" "," ",
                         paste("legend.muAD <- c(legend.muAD,expression(paste(mu[A*','*H],' = ',",
                               vector.of.muAD[k],",sep='')))",sep="")
  )))
}

B.star <- matrix(NaN,length(A_bar),length(vector.of.muAD))
E      <- matrix(NaN,length(A_bar),length(vector.of.muAD))
ExpEqH <- matrix(NaN,length(A_bar),length(vector.of.muAD))
PD.P   <- matrix(NaN,length(A_bar),length(vector.of.muAD))
PD.Q   <- matrix(NaN,length(A_bar),length(vector.of.muAD))
return.Eq.yearly <- matrix(NaN,length(A_bar),length(vector.of.muAD))

indic.plot <- 0

for(muAD in vector.of.muAD){
  
  indic.plot <- indic.plot + 1
  
  #mu_A$muprice_1[indic_D] <- muAD
  mu_A$muprice_1[indic_H] <- muAD
  
  # Update model accordingly:
  new_model_sol <- update.model_sol.4.mu_altern(model_sol,
                                                mu_altern = mu_A,
                                                indic_Euler = 1,
                                                H_if_Euler = horizon)
  
  #adjusted_A_bar <- A_bar * RF.strategy[horizon]
  adjusted_A_bar <- A_bar
  
  # Price risky bond:
  xx.fast <- varphi.hat.fast(new_model_sol,omega = omega.ZC,
                             H = horizon,x = x,
                             a = - omega_A,b = - log(adjusted_A_bar))
  xx2.fast <- varphi.bar.fast(new_model_sol,omega = omega_A,
                              H = horizon,x = x,
                              a = omega_A,b = log(adjusted_A_bar))
  B.star[,indic.plot] <- xx.fast + 1/adjusted_A_bar * xx2.fast
  
  # Price equity:
  xx.fast.E <- varphi.hat.fast(new_model_sol,omega = omega_A,
                               H = horizon,x = x,
                               a = - omega_A,b = - log(adjusted_A_bar))
  E[,indic.plot] <- xx.fast.E - adjusted_A_bar*xx.fast
  
  # Expected value of equity at maturity:
  EA1 <- fourier_complete(new_model_sol,x,omega_A,-omega_A,-log(adjusted_A_bar),h=horizon)
  EA2 <- fourier_complete(new_model_sol,x,0*omega_A,-omega_A,-log(adjusted_A_bar),h=horizon)
  ExpEqH[,indic.plot] <- EA1 - adjusted_A_bar*EA2
  
  # log Expected return on Equity:
  return.Eq.yearly[,indic.plot] <- log(ExpEqH[,indic.plot]/E[,indic.plot])/(model$tstep*horizon)
  
  # Default probabilities under P and Q:
  PD.P[,indic.plot] <- 1 - EA2
  PD.Q[,indic.plot] <- 1 - xx.fast/prices.ZCRF.bonds$P.t[horizon]
  
}



#Plots
FILE = paste("/outputs/Figures/Figure_Merton2.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=7, height=6)

par(mfrow=c(2,2))

par(plt=c(.15,.95,.2,.85))

xlab <- expression(bar(A))


# Plot of expected equity return
for(indic.plot in 1:length(vector.of.muAD)){
  if(indic.plot==1){
    plot(A_bar,100*return.Eq.yearly[,indic.plot],type="l",
         xlim=xlim,ylim=ylim.returns,xlab=xlab,
         ylab="",
         lwd=2,cex.lab=1.3,
         main="(a) Expected equity excess return (%)",
         col="black",las=1)
    legend("topleft",
           legend=legend.muAD,
           lty=c(1,3),
           col=c("black"),cex=1.1,
           lwd=c(2,2),bty = "n")
    # legend("bottomleft",
    #        legend=c("Indebted firm","Non-indedbted firm"),
    #        lty=c(1,1),
    #        col=c(P.col.line,"grey"),cex=1.1,
    #        lwd=c(2,2),bty = "n")
    #abline(h=100*return.Eq.yearly[1,indic.plot],lty=1,col="grey",lwd=2)
  }else{
    lines(A_bar,100*return.Eq.yearly[,indic.plot],
          lwd=2,lty=3,col="black")
    #abline(h=100*return.Eq.yearly[1,indic.plot],lty=3,col="grey",lwd=2)
  }
}

# Plot of yields
for(indic.plot in 1:length(vector.of.muAD)){
  if(indic.plot==1){
    plot(A_bar,-100/(horizon*model_sol$tstep)*log(B.star[,1]),
         lwd=2,type="l",
         xlim=xlim,ylim=ylim.ratioRates,xlab=xlab,
         ylab="",las=1,
         main="(b) Bond yields",
         cex.lab=1.3)
    legend("topleft",
           legend=c(legend.muAD,expression("risk-free yield")),
           lty=c(1,3,1),
           col=c("black","black","grey"),cex=1.1,
           lwd=c(2,2,2),bty = "n")
  }else{
    lines(A_bar,-100/(horizon*model_sol$tstep)*log(B.star[,2]),
          lwd=2,lty=3)
  }
}
abline(h=prices.ZCRF.bonds$r.t[horizon],col="grey",lwd=2)

# Plot of PDs
for(indic.plot in 1:length(vector.of.muAD)){
  if(indic.plot==1){
    plot(A_bar,log(PD.P[,indic.plot]),
         xlim=xlim,ylim=ylim.PD,lwd=2,type="l",
         xlab=xlab,yaxt="n",
         ylab="",
         main="(c) Probability of default (PD)",
         cex.lab=1.3,col="white")
    values.proba.axis <- c(.001,.005,.01,.05,.1,.25,.5,1)
    axis(side=2,
         at=log(values.proba.axis),
         labels=paste(100*values.proba.axis,"%",sep=""),las=1)
    for(i in 1:length(values.proba.axis)){
      abline(h=log(values.proba.axis[i]),lty=3,col="grey")
    }
    lines(A_bar,log(PD.P[,indic.plot]),col=P.col.line,lwd=2)
    lines(A_bar,log(PD.Q[,indic.plot]),col=Q.col.line,lwd=2)
    abline(h=0,col="grey")
    legend("topleft",
           legend=legend.muAD,
           lty=c(1,3),
           col=c("black"),cex=1.1,
           lwd=c(2,2),bty = "n")
    legend("bottomright",
           legend=c("Physical (P)","Risk-neutral (Q)"),
           lty=c(1,1),
           col=c(P.col.line,Q.col.line),cex=1.1,
           lwd=2,bty = "n")
  }else{
    lines(A_bar,log(PD.P[,indic.plot]),lwd=2,lty=3,col=P.col.line)
    lines(A_bar,log(PD.Q[,indic.plot]),col=Q.col.line,lwd=2,lty=3)
  }
}

# Plot of P/Q PD ratio
for(indic.plot in 1:length(vector.of.muAD)){
  if(indic.plot==1){
    plot(A_bar,PD.Q[,indic.plot]/PD.P[,indic.plot],
         lwd=2,type="l",
         xlim=xlim,ylim=ylim.ratioPD,xlab=xlab,
         ylab="",las=1,
         main="(d) Ratio of PDs (Q/P)",
         cex.lab=1.3)
    legend("topleft",
           legend=legend.muAD,
           lty=c(1,3),
           col=c("black"),cex=1.1,
           lwd=c(2,2),bty = "n")
  }else{
    lines(A_bar,PD.Q[,indic.plot]/PD.P[,indic.plot],lwd=2,lty=3)
  }
}



dev.off()


