
#Climate beta
y.lim   <- c(-.3,1.2)
y.limscc<- c(0,4000)
y.limzcb<- c(-5,5)
cexaxs  <- 1
cexlab  <- 1
cexmain <- 1

# Maturity, in periods, of the ZCB:
zcbmat  <-10

# Prepare omega vectors (for pricing):
omega_ZCB <- matrix(0,model_sol$n.X)
omega_T.at <- omega_ZCB
omega_T.at[which(model_sol$names.var.X=="T_at")] <- 1

# Considered values of mu_N (select three values):
values.of.mu_N <- c(0.001,round(model_sol$parameters$mu_N),
                    1.5*round(model_sol$parameters$mu_N))

# Considered values of mu_D:
seq.mu_D<-seq(0.001,
              2*model_sol$parameters$mu_D,
              #0.15,
              len=48)

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
                     
                     model_solswp1<-model_solve(modelswp1,theta0)
                     model_solswp2<-model_solve(modelswp2,theta0)
                     model_solswp3<-model_solve(modelswp3,theta0)
                     
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


splineswaps50 <-splinefun(seq.mu_D[!is.na(swaps50)],swaps50[!is.na(swaps50)],method = "hyman")
splineswaps0  <-splinefun(seq.mu_D[!is.na(swaps0)],swaps0[!is.na(swaps0)],method = "hyman")
splineswaps   <-splinefun(seq.mu_D[!is.na(swaps)],swaps[!is.na(swaps)],method = "hyman")

splinescc50   <-splinefun(seq.mu_D[!is.na(scc50)],scc50[!is.na(scc50)],method = "hyman")
splinescc0    <-splinefun(seq.mu_D[!is.na(scc0)],scc0[!is.na(scc0)],method = "hyman")
splinescc     <-splinefun(seq.mu_D[!is.na(scc)],scc[!is.na(scc)],method = "hyman")

splinezcb50   <-splinefun(seq.mu_D[!is.na(zcb50)],zcb50[!is.na(zcb50)],method = "hyman")
splinezcb0    <-splinefun(seq.mu_D[!is.na(zcb0)],zcb0[!is.na(zcb0)],method = "hyman")
splinezcb     <-splinefun(seq.mu_D[!is.na(zcb)],zcb[!is.na(zcb)],method = "hyman")

grid          <-seq(seq.mu_D[1],tail(seq.mu_D,1),length.out = 200+1)
grid          <-seq(seq.mu_D[1],tail(seq.mu_D,1),length.out = length(seq.mu_D))

swaps.50<-splineswaps50(grid)
swaps.0 <-splineswaps0(grid)
swaps.n <-splineswaps(grid)

scc.50  <-splinescc50(grid)
scc.0   <-splinescc0(grid)
scc.n   <-splinescc(grid)

zcb.50  <-splinezcb50(grid)
zcb.0   <-splinezcb0(grid)
zcb.n   <-splinezcb(grid)

climate.beta<-cbind(grid,swaps.n,swaps.50,swaps.0,
                    scc.n,scc.50,scc.0,
                    zcb.n,zcb.50,zcb.0)

mun  <- round(model_sol$parameters$mu_N,1)
mun3 <- round(values.of.mu_N[3],1)


#Plots
FILE = paste("/outputs/Figures/Figure_Climate_Beta.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=8, height=6)
par(plt=c(.1,.95,.25,.85))
par(mfrow=c(2,1))

plot(grid,swaps.n,type="l",lwd=3,col="black",
     ylim=y.lim,cex.main=cexmain,cex.axis=cexaxs,
     cex.lab=cexlab,
     xlab=expression(mu[D]), main="(a) - Temperature Risk Premium, maturity: 2100",
     ylab=expression(paste("Swap Price minus E(", T[AT*','*2100],")",
                           sep="")),
     las=0)
abline(h=0,lty=3,col="black")
abline(v=model_sol$parameters$mu_D,lty=3,col="lightslategrey")
lines(grid,swaps.0,col="dark grey",
      lty=2,lwd=3)
lines(grid,swaps.50,col="dark grey",
      lty=4,lwd=3)
legend("topleft",legend=c(expression(paste(mu[N],"=0")),
                          bquote(mu[N]==.(mun) ~" (baseline)"),
                          bquote(mu[N]==.(mun3) ~"")),
       col=c("dark grey","black","dark grey"),lty=c(2,1,4),lwd=3,
       bty = "n",cex=1)

# plot(grid,scc.n,type="l",lwd=3,col="black",
#      ylim=y.limscc,cex.main=cexmain,cex.axis=cexaxs,
#      cex.lab=cexlab,
#      xlab=expression(mu[D]),main="(b) - Social Cost of Carbon",
#      ylab=expression(paste("SCC (in $/GTC)"),
#                      sep=""),
#      las=0)
# abline(h=0,lty=3,col="black")
# abline(v=model_sol$parameters$mu_D,lty=3,col="lightslategrey")
# lines(grid,scc.0,col="black",
#       lty=2,lwd=3)
# lines(grid,scc.50,col="black",
#       lty=4,lwd=3)

plot(grid,zcb.n,type="l",lwd=3,col="black",
     ylim=y.limzcb,cex.main=cexmain,cex.axis=cexaxs,
     cex.lab=cexlab,
     xlab=expression(mu[D]),main="(b) - Long-Term Real Interest Rate (50 years)",
     ylab=expression(paste("Real Rate (in percent)"),
                     sep=""),
     las=0)
abline(h=0,lty=3,col="black")
abline(v=model_sol$parameters$mu_D,lty=3,col="lightslategrey")
lines(grid,zcb.0,col="dark grey",
      lty=2,lwd=3)
lines(grid,zcb.50,col="dark grey",
      lty=4,lwd=3)
dev.off()

#Separated Plots
###
FILE = paste("/outputs/Figures/Figure_Climate_Beta_swaps.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=8,width=5, height=3)
par(mfrow=c(1,1))
plot(grid,swaps.n,type="l",lwd=2,col="black",
     ylim=y.lim,
     xlab=expression(mu[D]), main="Temperature Risk Premium, maturity: 2100",
     ylab=expression(paste("Swap Price minus Expect. ", T[AT*','*2100]),
                     sep=""),
     las=0)
abline(h=0,lty=3,col="black")
abline(v=model_sol$parameters$mu_D,lty=3,col="lightslategrey")
lines(grid,swaps.0,col="dark grey",
      lty=2,lwd=2)
lines(grid,swaps.50,col="dark grey",
      lty=4,lwd=2)
legend("topleft",legend=c(expression(paste(mu[N],"=0")),
                          bquote(mu[N]==.(mun) ~" (baseline)"),
                          bquote(mu[N]==.(mun3) ~"")),
       col=c("dark grey","black","dark grey"),lty=c(2,1,4),lwd=rep(2),
       bty = "n",cex=1)
dev.off()
##
FILE = paste("/outputs/Figures/Figure_Climate_Beta_SCC.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=8,width=5, height=3)
par(mfrow=c(1,1))
plot(grid,scc.n,type="l",lwd=2,col="black",
     ylim=y.limscc,
     xlab=expression(mu[D]),main="Social Cost of Carbon",
     ylab=expression(paste("SCC (in $/GTC)"),
                     sep=""),
     las=0)
abline(h=0,lty=3,col="black")
abline(v=model_sol$parameters$mu_D,lty=3,col="lightslategrey")
lines(grid,scc.0,col="dark grey",
      lty=2,lwd=2)
lines(grid,scc.50,col="dark grey",
      lty=4,lwd=2)
legend("topleft",legend=c(expression(paste(mu[N],"=0")),
                          bquote(mu[N]==.(mun) ~" (baseline)"),
                          bquote(mu[N]==.(mun3) ~"")),
       col=c("dark grey","black","dark grey"),lty=c(2,1,4),lwd=rep(2),
       bty = "n",cex=1)
dev.off()
#
FILE = paste("/outputs/Figures/Figure_Climate_Beta_ZCB.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=8,width=5, height=3)
par(mfrow=c(1,1))
plot(grid,zcb.n,type="l",lwd=2,col="black",
     ylim=y.limzcb,
     xlab=expression(mu[D]),main="Long-term Real Interest Rate (50 years)",
     ylab=expression(paste("Real Rate (in percent)"),
                     sep=""),
     las=0)
abline(h=0,lty=3,col="black")
abline(v=model_sol$parameters$mu_D,lty=3,col="lightslategrey")
lines(grid,zcb.0,col="dark grey",
      lty=2,lwd=2)
lines(grid,zcb.50,col="dark grey",
      lty=4,lwd=2)
legend("topleft",legend=c(expression(paste(mu[N],"=0")),
                          bquote(mu[N]==.(mun) ~" (baseline)"),
                          bquote(mu[N]==.(mun3) ~"")),
       col=c("dark grey","black","dark grey"),lty=c(2,1,4),lwd=rep(2),
       bty = "n",cex=1)
dev.off()