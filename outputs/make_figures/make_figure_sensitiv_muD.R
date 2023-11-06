# ==============================================================================
# Figure illustrating sensitivity of temperature premium and long-term
# rates to mu_D and mu_N.
# ==============================================================================

y.lim   <- c(-.3,3.5)
y.limzcb<- c(-5,5)
cexaxs  <- 1
cexlab  <- 1
cexmain <- 1

len <- 32 # number of values of mu_D

# Maturity, in periods, of the ZCB:
zcbmat  <-10

# Prepare omega vectors (for pricing):
omega_ZCB <- matrix(0,model_sol$n.X)
omega_T.at <- omega_ZCB
omega_T.at[which(model_sol$names.var.X=="T_at")] <- 1

# Considered values of mu_N (select three values):
values.of.mu_N <- c(0.001,round(model_sol$parameters$mu_N),
                    2*round(model_sol$parameters$mu_N))

# Considered values of mu_D:
min.mu.D <- 0.001
max.mu.D <- 4*model_sol$parameters$mu_D
alpha <- .2

seq.mu_D <- (2*min.mu.D - max.mu.D) +
  2*(max.mu.D - min.mu.D) * exp(alpha*(0:(len-1)))/(1+exp(alpha*(0:(len-1))))
#plot(seq.mu_D)

# Run parallel computations:
cl <- makeCluster(number.of.cores)
registerDoParallel(cl)
save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))
clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))
all_res <- foreach(i = 1:length(seq.mu_D), 
                   .combine=rbind) %dopar% {
                     
                     modelswp1<-model_sol
                     modelswp2<-model_sol
                     modelswp3<-model_sol
                     
                     modelswp1$parameters$mu_D<-seq.mu_D[i]
                     modelswp2$parameters$mu_D<-seq.mu_D[i]
                     modelswp3$parameters$mu_D<-seq.mu_D[i]
                     
                     modelswp1$parameters$mu_N <- values.of.mu_N[1]
                     modelswp2$parameters$mu_N <- values.of.mu_N[2]
                     modelswp3$parameters$mu_N <- values.of.mu_N[3]
                     
                     model_solswp1<-model_solve(modelswp1)
                     model_solswp2<-model_solve(modelswp2)
                     model_solswp3<-model_solve(modelswp3)
                     
                     # Get values for model1:
                     swaps1    <-(varphi.tilde(model_solswp1,omega_T.at,H)[[1]]/
                                    varphi(model_solswp1,omega_ZCB,H)[[3]])[H]-
                       EV.fct(model_solswp1,H)$EX$T_at[H]
                     scc1      <- scc.fct(model_solswp1,0)
                     zcb1      <-varphi(model_solswp1,omega_ZCB,zcbmat)[[4]][zcbmat]                                          
                     # Get values for model2:
                     swaps2    <-(varphi.tilde(model_solswp2,omega_T.at,H)[[1]]/
                                    varphi(model_solswp2,omega_ZCB,H)[[3]])[H]-
                       EV.fct(model_solswp2,H)$EX$T_at[H]
                     scc2      <- scc.fct(model_solswp2,0)
                     zcb2      <-varphi(model_solswp2,omega_ZCB,zcbmat)[[4]][zcbmat]                                          
                     # Get values for model3:
                     swaps3    <-(varphi.tilde(model_solswp3,omega_T.at,H)[[1]]/
                                    varphi(model_solswp3,omega_ZCB,H)[[3]])[H]-
                       EV.fct(model_solswp3,H)$EX$T_at[H]
                     scc3      <- scc.fct(model_solswp3,0)
                     zcb3      <-varphi(model_solswp3,omega_ZCB,zcbmat)[[4]][zcbmat]                                          
                     
                     c(swaps1,swaps2,swaps3,scc1,scc2,scc3,zcb1,zcb2,zcb3)
                   }

stopCluster(cl)
file.remove("outputs/toto.Rdata")


all_res[abs(all_res)>10^5] <- NaN

swaps     <-all_res[,2]
scc       <-all_res[,5]
zcb       <-all_res[,8]

swaps0    <-all_res[,1]
scc0      <-all_res[,4]
zcb0      <-all_res[,7]

swaps50   <-all_res[,3]
scc50     <-all_res[,6]
zcb50     <-all_res[,9]


mun  <- round(model_sol$parameters$mu_N,1)
mun3 <- round(values.of.mu_N[3],1)

#Plots
FILE = paste("/outputs/Figures/Figure_sensitiv_muT.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=8, height=6)
par(plt=c(.1,.95,.25,.85))
par(mfrow=c(2,1))

plot(seq.mu_D,swaps,type="l",lwd=3,col="black",
     ylim=y.lim,cex.main=cexmain,cex.axis=cexaxs,
     cex.lab=cexlab,
     xlab=expression(paste("Damage uncertainty, ",mu[D],sep="")),
     main="(a) - Temperature Risk Premium, maturity: 2100",
     ylab=expression(paste("Swap Price minus E(", T[AT*','*2100],")",
                           sep="")),
     las=0)
abline(h=0,lty=3,col="black")
abline(v=model_sol$parameters$mu_D,lty=3,col="lightslategrey")
lines(seq.mu_D,swaps0,col="dark grey",
      lty=2,lwd=3)
lines(seq.mu_D,swaps50,col="dark grey",
      lty=4,lwd=3)
legend("topleft",legend=c(expression(paste(mu[N],"=0")),
                          bquote(mu[N]==.(mun) ~" (baseline)"),
                          bquote(mu[N]==.(mun3) ~"")),
       col=c("dark grey","black","dark grey"),lty=c(2,1,4),lwd=3,
       bty = "n",cex=1)

plot(seq.mu_D,zcb,type="l",lwd=3,col="black",
     ylim=y.limzcb,cex.main=cexmain,cex.axis=cexaxs,
     cex.lab=cexlab,
     xlab=expression(paste("Damage uncertainty, ",mu[D],sep="")),
     main="(b) - Long-Term Real Interest Rate (50 years)",
     ylab=expression(paste("Real Rate (in percent)"),
                     sep=""),
     las=0)
abline(h=0,lty=3,col="black")
abline(v=model_sol$parameters$mu_D,lty=3,col="lightslategrey")
lines(seq.mu_D,zcb0,col="dark grey",
      lty=2,lwd=3)
lines(seq.mu_D,zcb50,col="dark grey",
      lty=4,lwd=3)
dev.off()

