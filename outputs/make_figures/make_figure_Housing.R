# ==============================================================================
# Figure illustrating housing prices
# ==============================================================================

# Prepare omega vectors (for pricing):
omega_ZCB <- matrix(0,model_sol$n.X)
omega_T.at <- omega_ZCB
omega_T.at[which(model_sol$names.var.X=="T_at")] <- 1

# For Fourier transform:
x <- exp(seq(-5,5,length.out = 1000)) #grid for Proposition 8 (Fourier)

values.b_sk    <- c(model$parameters$b_sk,2*model$parameters$b_sk)
values.of.Hbar <- seq(0.4,4.0,length.out=40)


# Maturity/horizon:
H <- 98
#H <- 20

# Selection vector for R:
omega_R <- matrix(0,model_sol$n.X,1)
omega_R[which(model_sol$names.var.X=="Cum_dc")] <- 1

# Selection vector for H:
omega_H <- matrix(0,model_sol$n.X,1)
omega_H[which(model_sol$names.var.X=="H")] <- 1


# R_t's specification:
# mu_R <- list(muprice_0 = model_sol$target_vector[["mu_c0"]],
#              muprice_1 = matrix(0,model_sol$n.X,1))
mu_R <- list(muprice_0 = .02,
             muprice_1 = matrix(0,model_sol$n.X,1))


legend.b_sk <- NULL
for(k in 1:length(values.b_sk)){
  if(k==1){
    eval(parse(text = gsub(" "," ",
                           paste("legend.b_sk <- c(legend.b_sk,expression(paste(b[sk],' = ',",
                                 toString(values.b_sk[k]),",' (baseline)',sep='')))",sep="")
    )))
  }else{
    eval(parse(text = gsub(" "," ",
                           paste("legend.b_sk <- c(legend.b_sk,expression(paste(b[sk],' = ',",
                                 toString(values.b_sk[k]),",sep='')))",sep="")
    )))
  }
}
legend.b_sk[1]

all.house.prices.P <- matrix(NaN,length(values.of.Hbar),length(values.b_sk))
all.house.prices.Q <- matrix(NaN,length(values.of.Hbar),length(values.b_sk))

cl <- makeCluster(number.of.cores)
registerDoParallel(cl)

save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))

clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))

indic.values.bsk <- 0
for(b_sk in values.b_sk){
  indic.values.bsk <- indic.values.bsk + 1
  
  print(paste("  Value of b_sk: ",indic.values.bsk," out of ",length(values.b_sk),sep=""))
  
  if(b_sk==model$parameters$b_sk){
    new_model_sol <- model_sol
  }else{
    new_model <- model
    new_model$parameters$b_sk <- b_sk
    new_model_sol<-model_solve(new_model)
  }
  
  # Returns of risk-free strategy:
  Price.ZC <- varphi(new_model_sol,omega_ZCB,H = H)[["P.t"]]
  
  new_model_sol <- update.model_sol.4.mu_altern(new_model_sol,
                                                mu_altern = mu_R,
                                                indic_Euler = 0,
                                                H_if_Euler = NaN)
  
  # Compute proba that H_t > Hbar
  
  all.Probas.P <- foreach(h = 1:H, .combine=cbind) %dopar% {
    probas <- fourier(new_model_sol,x,values.of.Hbar,h,
                      which(model_sol$names.var.X=="H"))
    probas
  }
  print("   *   Physical proba: done")
  all.discounted.payoffs.Q <- foreach(h = 1:H, .combine=cbind) %dopar% {
    NPV<-c(varphi.hat.fast(new_model_sol,omega = omega_R,
                           H=h,x,a = omega_H,
                           b=values.of.Hbar))
    NPV
  }
  print("   **  NPV: done")
  
  # House prices:
  house.prices.Q <- apply(all.discounted.payoffs.Q,1,sum)
  
  # House prices under P:
  expected.R <- apply(matrix(1:H,ncol=1),1,
                      function(h){multi.lt.fct.N(new_model_sol,U=omega_R,h=h)})
  discounted.R <- matrix(1,length(values.of.Hbar),1) %*% 
    matrix(expected.R*Price.ZC,nrow=1)
  expected.payoffs <- discounted.R * all.Probas.P
  house.prices.P <- apply(expected.payoffs,1,sum)
  
  all.house.prices.P[,indic.values.bsk] <- house.prices.P * model_sol$tstep
  all.house.prices.Q[,indic.values.bsk] <- house.prices.Q * model_sol$tstep
}


stopCluster(cl)
file.remove("outputs/toto.Rdata")

par(mfrow=c(2,2))
plot(discounted.R[1,],type="l")
plot(all.Probas.P[10,],type="l")


FILE = paste("/outputs/Figures/Figure_Housing.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=11,width=9, height=4)

par(mfrow=c(1,2))
par(plt=c(.16,.95,.2,.85))

plot(values.of.Hbar,all.house.prices.Q[,1],col=Q.col.line,type="l",lwd=2,
     main="(a) House prices",xlab=expression(bar(H)),
     ylab="in multiples of initial year's rental income",
     ylim=c(min(all.house.prices.Q,all.house.prices.P),
            max(all.house.prices.Q,all.house.prices.P)),
     las=1)
lines(values.of.Hbar,all.house.prices.P[,1],col=P.col.line,lwd=2,lty=1)
lines(values.of.Hbar,all.house.prices.Q[,2],col=Q.col.line,lwd=2,lty=2)
lines(values.of.Hbar,all.house.prices.P[,2],col=P.col.line,lwd=2,lty=2)
legend("bottomright",
       legend=c("House price (incl. risk premiums)","Discounted expected payoffs"),
       lty=c(1,1),
       col=c(Q.col.line,P.col.line),cex=1,
       lwd=c(2,2),bty = "n")

legend("topleft",
       legend=legend.b_sk,
       lty=c(1,2),
       col="black",cex=1,
       lwd=c(2,2),bty = "n")


shares.RP <- 100*(all.house.prices.P/all.house.prices.Q-1)

plot(values.of.Hbar,shares.RP[,1],type="l",lwd=2,
     ylim=c(0,max(shares.RP)),
     main="(b) Share of risk premiums",xlab=expression(bar(H)),ylab="in percent",las=1)
lines(values.of.Hbar,shares.RP[,2],lwd=2,lty=2)

dev.off()
