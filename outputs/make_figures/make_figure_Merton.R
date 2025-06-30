# ==============================================================================
# FIGURE VI.3. Merton model
# Figure_Merton2.pdf
# Figure_Merton1.pdf (not reported in Supplemental Appendix)
# ==============================================================================

# Figures illustrating Merton model (1/2)
ALPHA <- .0

nb.values.variable <- 2000 # for extrapolation

y.lim <- c(.5,14)
x.lim <- c(2030,2070)

# For Fourier transform:
x <- exp(seq(-10,10,length.out = 1000)) #grid for Proposition 8 (Fourier)

values.of.logA <- log(seq(.01,200,length.out=1000))


# Maturity/horizon:
H <- (x.lim[2] - model_sol$vec_date[1])/model_sol$tstep

# Returns of risk-free strategy:
omega.ZC <- matrix(0,model_sol$n.X,1)
prices.ZCRF.bonds   <- varphi(model_sol,
                              omega.varphi = omega.ZC,
                              H = H)
RF.strategy <- 1/prices.ZCRF.bonds$P.t

# Selection vector for A:
omega_A <- matrix(0,model_sol$n.X,1)
omega_A[which(model_sol$names.var.X=="Cum_dc")] <- 1

vector.of.elasticity.wrt.dc <- c(2.5,2.5)
vector.of.muAD              <- c(0,-.2*model_sol$tstep)

# A_t's specification:
mu_A <- list(muprice_0 = 0,
             muprice_1 = matrix(0,model_sol$n.X,1))

indic_delc <- which(model_sol$names.var.X=="delc")

# in case one wants to add idiosyncratic shocks:
indic_X <- which(model_sol$names.var.X=="eta_X")
mu_A$muprice_1[indic_X] <- .0

indic_D    <- which(model_sol$names.var.X=="D")
indic_H    <- which(model_sol$names.var.X=="H")
indic_delc <- which(model_sol$names.var.X=="delc")

panel.titles <- NULL
for(k in 1:length(vector.of.muAD)){
  eval(parse(text = gsub(" "," ",
                         paste("panel.titles <- c(panel.titles,expression(paste('Panel (",
                               letters[k],") ',mu[A*','*D],' = ',",
                               vector.of.muAD[k],",sep='')))",sep="")
  )))
}


all.expected.A <- NULL

# ------------------------------------------------------------------------------
# Plots----
# 1----
FILE = paste("/outputs/Figures/Figure_Merton1.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=8, height=4)

par(mfrow=c(1,length(vector.of.muAD)))
par(plt=c(.15,.95,.15,.85))

indic.plot <- 0

cl <- makeCluster(number.of.cores)
registerDoParallel(cl)

save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))

clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))

for(muAD in vector.of.muAD){
  
  indic.plot <- indic.plot + 1
  
  iiii <- which(muAD==vector.of.muAD)
  mu_A <- list(muprice_0 = 0,
               muprice_1 = matrix(0,model_sol$n.X,1))
  mu_A$muprice_1[indic_delc] <- vector.of.elasticity.wrt.dc[iiii]
  mu_A$muprice_1[indic_H] <- muAD
  if(muAD==0){
    mu_A$muprice_1[indic_X] <- ALPHA
  }
  
  # Update model accordingly:
  model_sol <- update.model_sol.4.mu_altern(model_sol,
                                            mu_altern = mu_A,
                                            indic_Euler = 1,
                                            H_if_Euler = H)
  
  expected.A <- apply(matrix(1:H,ncol=1),1,
                      function(h){multi.lt.fct.N(model_sol,U=omega_A,h=h)})
  
  # save results (for comparison):
  all.expected.A <- rbind(all.expected.A,
                          expected.A)
  
  # Compute P and Q proba:
  Price.ZC <- varphi(model_sol,omega.ZC,H = H)[[3]]
  all.Probas.P <- foreach(h = 1:H, .combine=cbind) %dopar% {
    probas <- fourier(model_sol,x,values.of.logA,h,
                      which(model_sol$names.var.X=="Cum_dc"))
    probas
  }
  
  # Compute confidence intervals:
  CI.P <- confidence_intervals_across_horizons(all.Probas.P,
                                               values.of.variable = exp(values.of.logA),
                                               nb.values.variable = nb.values.variable,
                                               vector.of.CI = vector.of.CI)
  
  scale.A.values <- CI.P$scale.variable.values
  all.CI.P  <- CI.P$all.CI
  all.pdf.P <- CI.P$all.pdf
  all.cdf.P <- CI.P$all.cdf
  
  plot(model_sol$vec_date[2:(H+1)],RF.strategy,
       xlim=x.lim,ylim=y.lim, cex.main=1.5,cex.axis=1.5,cex.lab=1.5,
       col="white",xlab="",ylab="",las=1,
       main=panel.titles[indic.plot])
  
  for(i in length(vector.of.CI):1){
    # P
    polygon(c(model_sol$vec_date[2:(H+1)],rev(model_sol$vec_date[2:(H+1)])),
            c(all.CI.P[1,,i],rev(all.CI.P[2,,i])),
            col=P.col,border = NA)
  }
  lines(model_sol$vec_date[2:(H+1)],
        expected.A,lwd=2,col=P.col.line)
  lines(model_sol$vec_date[2:(H+1)],RF.strategy,
        lty=2,lwd=2,col="lightsteelblue4")
  
  grid()
  
  if(indic.plot==1){
    legend("topleft",
           legend=c("Expected value of firm's assets","Risk-free return"),
           lty=c(1,2),
           col=c(P.col.line,"lightsteelblue4"),cex=1.4,
           lwd=c(2,2),bty = "n")
  }
}

dev.off()

stopCluster(cl)
file.remove("outputs/toto.Rdata")

#print(log(all.expected.A[,10])/50)



ylim.returns <- 100*c(.02,.045)
ylim.PD      <- c(log(.0001),log(.25))
ylim.ratioPD <- c(0,10)
ylim.ratioRates <- c(2,3.5)
xlim         <- c(.1,1.2)

# Values of debt to repay:
A_bar <- seq(.01,max(1,max(xlim)),length.out=50)

horizon <- 4 # in number of periods

indic_D    <- which(model_sol$names.var.X=="D")
indic_H    <- which(model_sol$names.var.X=="H")
indic_delc <- which(model_sol$names.var.X=="delc")

# For Fourier-transform computation:
x <- exp(seq(-5,5,length.out = 10000)) # grid for Fourier

# Returns of risk-free strategy:
omega.ZC <- matrix(0,model_sol$n.X,1)
prices.ZCRF.bonds   <- varphi(model_sol,
                              omega.varphi = omega.ZC,
                              H = horizon)
RF.strategy <- 1/prices.ZCRF.bonds$P.t

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
  
  iiii <- which(muAD==vector.of.muAD)
  mu_A <- list(muprice_0 = 0,
               muprice_1 = matrix(0,model_sol$n.X,1))
  mu_A$muprice_1[indic_delc] <- vector.of.elasticity.wrt.dc[iiii]
  mu_A$muprice_1[indic_H] <- muAD
  if(muAD==0){
    mu_A$muprice_1[indic_X] <- ALPHA
  }
  
  # Update model accordingly:
  new_model_sol <- update.model_sol.4.mu_altern(model_sol,
                                                mu_altern = mu_A,
                                                indic_Euler = 1,
                                                H_if_Euler = horizon)
  
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
  return.Eq.yearly[,indic.plot] <- log(ExpEqH[,indic.plot]/E[,indic.plot])/(model_sol$tstep*horizon)
  
  # Default probabilities under P and Q:
  PD.P[,indic.plot] <- 1 - EA2
  PD.Q[,indic.plot] <- 1 - xx.fast/prices.ZCRF.bonds$P.t[horizon]
  
}

# 2----
FILE = paste("/outputs/Figures/Figure_Merton2.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=8, height=4)

par(mfrow=c(1,2))

par(plt=c(.15,.95,.2,.85))

xlab <- expression(bar(A))

# Plot of yields
for(indic.plot in 1:length(vector.of.muAD)){
  if(indic.plot==1){
    plot(A_bar,-100/(horizon*model_sol$tstep)*log(B.star[,1]),
         lwd=2,type="l",
         xlim=xlim,ylim=ylim.ratioRates,xlab=xlab,
         ylab="",las=1,
         main="(a) Bond yields",
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
         main="(b) Probability of default (PD)",
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

dev.off()





