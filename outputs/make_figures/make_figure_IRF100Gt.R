# ==============================================================================
# Response of temperature to a 1-GtC shock (M_at)
# ==============================================================================

indic.M_at <- which(model_sol$names.var.X=="M_at")

# Load EPA data:
load(file="/Users/jrenne/Dropbox/Research/TIBs/CR_Rcodes/data/Figure223_EPA_data.Rdat")

shock <- 1 # Gt Carbon

model_sol_shock <- model_sol
model_sol_shock$vector.ini$ini_Mat <- model_sol_shock$vector.ini$ini_Mat + shock
model_sol_shock$X[indic.M_at]      <- model_sol_shock$vector.ini$ini_Mat

end.date <- 2300
h.end <- (end.date-model_sol$vec_date[1])/model_sol$tstep - 1
EV       <- EV.fct(model_sol,h=h.end)
EV_shock <- EV.fct(model_sol_shock,h=h.end)

IRF <- EV_shock$EX$M_at - EV$EX$M_at
#plot(IRF,type="l")
IRF <- EV_shock$EX$T_at - EV$EX$T_at

lower.bound.01 <- fair$q01
upper.bound.99 <- fair$q99
lower.bound.05 <- fair$q05
upper.bound.95 <- fair$q95


FILE = "/outputs/Figures/Figure_IRF1GtC.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=11, width=8, height=4)

par(mfrow=c(1,1))
par(plt=c(.15,.95,.15,.95))
#par(mar=c(5, 6, 4, 2))

plot(EV$date+5,IRF,col="white",las=1,
     xlim=c(2025,2300),
     ylim=c(0,1.1*max(upper.bound.99)),
     xlab="",
     ylab="")
grid()

title(ylab='Temperature Anomaly from 1GtC in 2030', line=5,
      cex.lab=1)

polygon(x=c(fair$year,rev(fair$year)),c(lower.bound.01,rev(upper.bound.99)),
        col='grey90',border=NaN)
polygon(x=c(fair$year,rev(fair$year)),c(lower.bound.05,rev(upper.bound.95)),
        col='grey80',border=NaN)

lines(EV$date+5,IRF,lwd=1,pch=3,type="b")
lines(fair$year,fair$mean,lwd=2,col="black")
lines(hector$year,hector$temp.delta,col="#E69F00",lty=3,lwd=3)
lines(magicc$year,magicc$temp.delta,col="#56B4E9",lty=2,lwd=3)

legend("topleft",
       legend=c("CR","HECTOR 2.5","MAGICC 7.5.3",
                           "FaIR 1.6.2 (mean)"),
       col=c("black","#E69F00", "#56B4E9","black"),
       lty=c(NaN,3,2,1,NaN,NaN),
       lwd=c(1,3,3,2),
       pch=c(3,NaN,NaN,NaN),
       bty = "n",cex=1,
       bg="white",
       ncol=2)

dev.off()
