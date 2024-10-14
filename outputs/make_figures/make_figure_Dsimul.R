# ==============================================================================
# Simulation of damages
# ==============================================================================

#Gamma 0 shock simulations
nb       <- 12
nb.simul <- model_sol$Tmax                                                                   

H <- model_sol$horiz.2100

simul_d <-simul.function(model_sol,nb.simul,nb)
simD <- simul_d$D[1:H,]
simN <- simul_d$N[1:H,]

EV <- EV.fct(model_sol,h=H)
ED <- EV$EX$D
EN <- EV$EX$N

#Plots
#D
FILE = "/outputs/Figures/Figure_Cum_Damages_simul.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=10, height=8)
par(mfrow=c(2,1))  
plot(model_sol$vec_date[2:(H+1)],ED,col=P.col.line,type="l",
     ylim=c(range(ED,simD)),
     lwd=2,lty=2,las=1,
     ylab="Consumption Damages (in percent)", xlab="",
     cex.main=1.2,cex.axis=1.2,cex.lab=1.2)
for(i in 1:dim(simD)[2]){
  lines(model_sol$vec_date[2:(H+1)],simD[,i],
        col="grey75")
}
legend("bottomleft",legend=c("Expected Path","Simulations"),
       col=c(P.col.line,"grey75"),lty=c(2,1),lwd=c(2,1),
       bty = "n",cex=1.2)
plot(model_sol$vec_date[2:(H+1)],EN,col=P.col.line,type="l",
     ylim=c(range(EN,simN)),
     lwd=2,lty=2,las=1,
     ylab="Emissions from Permafrost", xlab="",
     cex.main=1.2,cex.axis=1.2,cex.lab=1.2)
for(i in 1:dim(simN)[2]){
  lines(model_sol$vec_date[2:(H+1)],simN[,i],
        col="grey75")
}
# lines(model_sol$vec_date[2:(H+1)],smallN,
#       col=P.col.line)
legend("topleft",legend=c("Expected Path","Simulations"),
       col=c(P.col.line,"grey75"),lty=c(2,1),lwd=c(2,1),
       bty = "n",cex=1.2)

dev.off()
