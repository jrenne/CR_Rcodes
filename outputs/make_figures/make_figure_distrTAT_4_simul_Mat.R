# ==============================================================================
# Figure 1*: From carbon concentrations to atmospheric temperature
# Figure_StochSimul_RCP_to_TAT.pdf
# ==============================================================================

#* Figure that illustrates the influence of the linearization of F versus T_AT
#* on temperatures.
#* The baseline model is used to simulate M_AT trajectories.
model_sol_aux <- model_sol

nb.t <- (2300-2020)/model_sol_aux$tstep # maximum horizon
nb.traj <- 40000 # number of simulated trajectories

# Step 1: Get simulated trajectories of M_AT (using baseline model):
res <- simul.function(model_sol_aux,nb.simul.t=nb.t,nb.traj=nb.traj)#,setseed=123)
EV <- EV.fct(model_sol_aux,h=nb.t)
Mat.trajectory <- res$M_at

# Step 2: Determine resulting T_AT paths, for both the linear and non-linear model:
res_CR    <- simul_TAT_condit_MAT(model_sol_aux,Mat.trajectory,indic_stochastic = TRUE)

deviations <- (res_CR$Tat.linear[nb.t,] - res_CR$Tat.nonlinear[nb.t,])/res_CR$Tat.nonlinear[nb.t,]
mean(abs(deviations))

dens_linear    <- density(res_CR$Tat.linear[nb.t,],bw = .1)
dens_nonlinear <- density(res_CR$Tat.nonlinear[nb.t,],bw = .1)

mean_linear    <- mean(res_CR$Tat.linear[nb.t,])
mean_nonlinear <- mean(res_CR$Tat.nonlinear[nb.t,])

stdv_linear    <- sd(res_CR$Tat.linear[nb.t,])
stdv_nonlinear <- sd(res_CR$Tat.nonlinear[nb.t,])

print(matrix(c(mean_linear,mean_nonlinear,stdv_linear,stdv_nonlinear),2,2))

ylim <- c(0,max(dens_nonlinear$y))

# ------------------------------------------------------------------------------
# Plot----
col_linear    <- Q.col.line
col_nonlinear <- P.col.line

FILE = "/outputs/Figures/Figure_StochSimul_RCP_to_TAT.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=6, height=4)

par(plt=c(.1,.95,.2,.95))

plot(dens_nonlinear,ylim=ylim,lwd=2,col=col_nonlinear,
     xlab="Atmospheric Temperature in 2300, in °C",ylab="",
     main="",las=1)
lines(dens_linear,col=col_linear,lwd=2)

abline(v=mean_linear,col=col_linear,lty=2)
abline(v=mean_nonlinear,col=col_nonlinear,lty=2)

points(mean_linear,.3*ylim[2],pch=19,col=col_linear,cex=1.5)
points(mean_nonlinear,.2*ylim[2],pch=19,col=col_nonlinear,cex=1.5)

lines(c(mean_linear - stdv_linear,
        mean_linear + stdv_linear),c(.3*ylim[2],.3*ylim[2]),
      col=col_linear,lwd=2,type="b",pch=3)
lines(c(mean_nonlinear - stdv_nonlinear,
        mean_nonlinear + stdv_nonlinear),c(.2*ylim[2],.2*ylim[2]),
      col=col_nonlinear,lwd=2,type="b",pch=3)

legend("topright",
       legend=c("Linearized specif. (= CR baseline)",
                "Non-linearized specif. (= CDICE)"),
       lty=c(1,1),
       col=c(col_linear,col_nonlinear),
       pch=c(NaN),
       lwd=c(2),
       seg.len = 3,
       bty = "n",cex=1)

legend("topleft",
       legend=c("Mean",
                "Mean +/- 1 std dev."),
       lty=NaN,
       pch=c(19,3),
       lwd=c(2),
       seg.len = 1,
       bty = "n",cex=1)

dev.off()