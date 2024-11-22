# ==============================================================================
# Figure illustrating Merton model (1/2)
# ==============================================================================

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
#vector.of.muAD              <- c(0,-elasticity.wrt.dc)
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

#Plots
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
  #mu_A$muprice_1[indic_D] <- muAD
  mu_A$muprice_1[indic_H] <- muAD
  if(muAD==0){
    mu_A$muprice_1[indic_X] <- ALPHA
  }
  
  # Update model accordingly:
  model_sol <- update.model_sol.4.mu_altern(model_sol,
                                            mu_altern = mu_A,
                                            indic_Euler = 1,
                                            H_if_Euler = H)
  
  # # Check Euler equations:
  # omega_A <- matrix(0,model_sol$n.X,1)
  # omega_A[which(model_sol$names.var.X=="Cum_dc")] <- 1 # where A_t is (tilde{A} at that stage).
  # prices.Atilde <- varphi(model_sol,omega.varphi = omega_A,H=H)
  # plot(prices.Atilde$P.t)
  
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


print(log(all.expected.A[,10])/50)

