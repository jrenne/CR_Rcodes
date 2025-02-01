# ==============================================================================
# Figure 2*: Simulated damages
# Figure_simul_cumD.pdf
# ==============================================================================

nb.t <- 20
nb.traj <- 10
nb.t.2.show <- 10

res <- simul.function(model_sol,nb.simul.t=nb.t,nb.traj=nb.traj)
EV <- EV.fct(model_sol,h=nb.t)

aux <- apply(res$D,2,cumsum)
CumDamages <- exp(aux) - 1
CumN <- apply(res$N,2,cumsum)


E_Cumdamages <- exp(cumsum(EV$EX$D)) - 1

E_CumN <- cumsum(EV$EX$N)

# ------------------------------------------------------------------------------
# Plot----
FILE = "/outputs/Figures/Figure_simul_cumD.pdf"
pdf(file=FILE,pointsize=12, width=11, height=5)

par(mfrow=c(1,2))
par(plt=c(.15,.95,.1,.85))

plot(EV$date,E_Cumdamages,col="red",ylim=c(0,max(CumDamages)),type="l",lwd=3,
     main="(a) Cumulated damages",
     ylab="Damages (fraction of pre-damage output)",xlab="",lty=2)
for(i in 1:nb.t.2.show){
  lines(EV$date,CumDamages[,i],col=i+1)
}



plot(EV$date,E_CumN,col="red",ylim=c(0,max(CumN)),type="l",lwd=3,
     main="(b) Cumulated permafrost-related carbon emissions",
     ylab="Permafrost-related carbon emissions (in GtCO2)",xlab="",lty=2)
for(i in 1:nb.t.2.show){
  lines(EV$date,CumN[,i],col=i+1)
}

dev.off()


