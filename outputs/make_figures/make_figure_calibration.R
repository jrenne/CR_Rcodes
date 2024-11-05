# ==============================================================================
# Figure illustrating the calibration approach
# ==============================================================================

# Plot ----
FILE = "/outputs/Figures/Figure_Calibration.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=11, width=7, height=6)
par(plt=c(.2,.95,.25,.85))

par(mfrow=c(2,2))

# Damages ----------------------------------------------------------------------
make_figure_calibration_Damages(model_sol,
                                main.title = "(a) 2100 Fractional damages")


# Sea level rise ---------------------------------------------------------------
make_figure_calibration_SLR(model_sol,
                            main.title = "(b) 2100 Sea-level rise")


# Permafrost-related emissions -------------------------------------------------
make_figure_calibration_N(model_sol,
                          CumN.values = seq(0,2000,length.out=1000),
                          nb.values.interpol=1000,
                          xN = exp(seq(-20,20,length.out = 5000)),
                          main.title="(c.1) 2100 Permafrost-related emissions")

all.Ttilde <- c(2,4,6)

CumN.values <- seq(0,3000,length.out=60)
CumN.values.interpol <- seq(CumN.values[1],
                            CumN.values[length(CumN.values)],
                            length.out=1000)
# to compute Riemann sum:
xN <- exp(seq(-20,20,length.out = 2000))

all.cdf <- matrix(NaN,length(CumN.values.interpol),length(all.Ttilde))
all.pdf <- matrix(NaN,length(CumN.values.interpol)-1,length(all.Ttilde))

for(i in 1:length(all.Ttilde)){
  aux <- Fourier.psi(model_sol,CumN.values,xN,
                     T0=model$vector.ini$ini_Tat,
                     Tstar = all.Ttilde[i],
                     tstar=NaN,
                     psi.LT=LT.CumN.infinite)
  
  tmp <- splinefun(x=CumN.values, y=aux, method="hyman")
  
  fitted.cdf.values <- tmp(CumN.values.interpol)
  fitted.pdf.values <- diff(fitted.cdf.values)
  
  all.cdf[,i] <- fitted.cdf.values
  all.pdf[,i] <- fitted.pdf.values
} 

xlim <- c(0,max(CumN.values))
plot(CumN.values.interpol[2:length(CumN.values.interpol)],
     all.pdf[,1],type="l",col="grey",las=1,yaxt="no",
     main="(c.2) Total permafrost-related emissions",
     xlab="GtC",ylab="",lwd=2,xlim=xlim,lty=2)
lines(CumN.values.interpol[2:length(CumN.values.interpol)],
      all.pdf[,2],lwd=2,lty=1)
lines(CumN.values.interpol[2:length(CumN.values.interpol)],
      all.pdf[,3],lwd=2,col="grey",lty=3)
abline(v=model_sol$target_vector[["ECumNinf"]],col="red",lwd=2)

legend("topright",
       legend=paste(all.Ttilde,"\u00B0C",sep=""),
       title=expression(paste("Long-term value of ",T[infinity],":",sep="")),
       lty=c(2,1,3),
       col=c("grey","black","grey"),
       #cex=1.5,
       lwd=2,bty = "n")

dev.off()
