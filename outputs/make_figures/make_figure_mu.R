x.lim <- c(model_sol$vec_date[2],2200)

dice_data           <- read.csv("data/mu.csv",sep=";",header = F)
row.names(dice_data)<- dice_data[,1]
dice_data           <- dice_data[-1]

mu_dice                    <-matrix(0,dim(dice_data)[2],1)
mu_dice[1:length(mu_dice)] <-apply(dice_data[2,],2,function(x)min(x,1))

#Plot
FILE = paste("/outputs/Figures/Figure_Mitigation_comparison.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=5, height=3)
plot(model_sol$vec_date[2:length(model_sol$vec_date)],
     model_sol$mu[1:(length(model_sol$vec_date)-1)],
     type="l", xlab="Year",ylab="",lwd=2,
     xlim=x.lim,las=1)
lines(model_sol$vec_date[2:length(model_sol$vec_date)],
      mu_dice[2:length(model_sol$vec_date)],lwd=2,lty=2,
      col="grey")
legend("topleft",
       title = "Mitigation rate:",
       legend=c("Present study","DICE"),
       lty=c(1,2),
       col=c("black","grey"),
       lwd=c(2,2),
       bty = "n",cex=1)
dev.off()
