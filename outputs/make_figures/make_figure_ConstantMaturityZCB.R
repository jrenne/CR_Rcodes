# ==============================================================================
# Figure showing expected trajectories of long-term rates
# ==============================================================================

EV <- EV.fct(model_sol)
Price.ZC <- varphi(model_sol,
                   omega.varphi=matrix(0,model_sol$n.X,1),
                   H = model_sol$horiz.2100)

FILE = paste("/outputs/Figures/Figure_cstH_ZCB.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=9, height=4)
par(plt=c(.1,.95,.2,.85))
par(mfrow=c(1,2))

maturities = c(2,10)

expected.yds <- compute_CMT(model_sol,EV,Price.ZC,maturities)

make_figure_CMT(expected.yds,model_sol,maturities)

dev.off()
