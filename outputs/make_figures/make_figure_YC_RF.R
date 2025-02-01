# ==============================================================================
# FIGURE 6. The term structure of real rates
# Figure_YC_RF.pdf
# FIGURE V.3. Comparison of temperature risk premiums
# Figure_TRP_comparison.pdf
# ==============================================================================

# Maximum maturity (in years):
H <- 300

col.Lemoine <- "#B22222"
col.BKO     <- "#008B8B"
col.CR      <- "black"
col.BR      <- "darkgoldenrod3"

# Lemoine (2021) model ---------------------------------------------------------

res <- make.Lemoine.model() # load baseline specification
model.Lemoine <- res$model
X0 <- res$X0

chge <- .4 # reflects uncertainty on damages and climate sensitivity

res <- compute.SCC.Lemoine(model.Lemoine,
                           damage.type = "G. Lemoine",
                           X0,H=H,nb.replic=200,
                           s.values = c(1-chge,1,1+chge)*model.Lemoine$s,
                           coef.multip.values = c(1-chge,1,1+chge),
                           seed0 = 123,
                           shock = 1)

yds_lemoine <- -log(apply(res$MM,1,mean))/1:H

# Temperature risk premium:
TRP_lemoine <- res$E.T.risk.adj - res$E.T


# BKO (2019) model -------------------------------------------------------------

model.BKO <- make.BKO.model()
model_sol.BKO <- solve_model.bko(model.BKO)
Tbar <- model.BKO$Tbar
# yields:
res_prices <- compute_ab.bko(model_sol.BKO,H)
yds_BKO <- res_prices$all_r_a + res_prices$all_r_b*Tbar
# Temperature risk premium:
epsilon <- .0001
B <- exp(res_prices$all_a + res_prices$all_b*Tbar)
res_prices_eps <- compute_ab.bko(model_sol.BKO,H,u=epsilon)
Tswap <- ((res_prices_eps$all_a - res_prices$all_a) +
            (res_prices_eps$all_b - res_prices$all_b)*Tbar)/epsilon
TRP_BKO <- Tswap - Tbar

# BKO with higher damages:
multip.factor <- 2
model.BKO$d <- multip.factor * model.BKO$d
model_sol.BKO <- solve_model.bko(model.BKO)
# yields:
res_prices <- compute_ab.bko(model_sol.BKO,H)
yds_BKO_highDamages <- res_prices$all_r_a + res_prices$all_r_b*Tbar

# Temperature risk premium:
B <- exp(res_prices$all_a + res_prices$all_b*Tbar)
res_prices_eps <- compute_ab.bko(model_sol.BKO,H,u=epsilon)
Tswap <- ((res_prices_eps$all_a - res_prices$all_a) +
            (res_prices_eps$all_b - res_prices$all_b)*Tbar)/epsilon
TRP_BKO_highDamages <- Tswap - Tbar


# CR model EZ ------------------------------------------------------------------

# Prepare omega vectors (for pricing):
omega_ZCB <- matrix(0,model_sol$n.X)
omega_T.at <- omega_ZCB
omega_T.at[which(model_sol$names.var.X=="T_at")] <- 1

Price.ZC <- varphi(model_sol,
                   omega.varphi=omega_ZCB,
                   H = H/model_sol$tstep)
yds_CR_EZ <- Price.ZC$r.t/100

EV   <- EV.fct(model_sol,h=H/model_sol$tstep)
ET.P_CR_EZ <- EV$EX$T_at[1:(H/model_sol$tstep)]
ET.Q_CR_EZ <- varphi.tilde(model_sol,omega_T.at,H/model_sol$tstep)[[1]]/
  varphi(model_sol,omega_ZCB,H/model_sol$tstep)[[3]]
TRP_CR_EZ <- ET.Q_CR_EZ - ET.P_CR_EZ

# CR model EZ (low delta) ------------------------------------------------------
model_lowDelta <- model_sol
aux <- - log(model_sol$parameters$delta)/model_sol$tstep
new_aux <- aux - .005
model_lowDelta$parameters$delta <- (1 - new_aux)^model_sol$tstep
model_lowDelta <- solveParam4c(model_lowDelta)
model_sol_lowDelta <- model_solve(model_lowDelta)

Price.ZC <- varphi(model_sol_lowDelta,
                   omega.varphi=omega_ZCB,
                   H = H/model_sol$tstep)
yds_CR_EZ_lowDelta <- Price.ZC$r.t/100

EV   <- EV.fct(model_sol_lowDelta,h=H/model_sol$tstep)
ET.P_CR_EZ <- EV$EX$T_at[1:(H/model_sol$tstep)]
ET.Q_CR_EZ <- varphi.tilde(model_sol_lowDelta,omega_T.at,H/model_sol$tstep)[[1]]/
  varphi(model_sol_lowDelta,omega_ZCB,H/model_sol$tstep)[[3]]
TRP_CR_EZ_lowDelta <- ET.Q_CR_EZ - ET.P_CR_EZ

# CR model EZ (high delta) ------------------------------------------------------
model_highDelta <- model_sol
aux <- - log(model_sol$parameters$delta)/model_sol$tstep
new_aux <- aux + .005
model_highDelta$parameters$delta <- (1 - new_aux)^model_sol$tstep
model_highDelta <- solveParam4c(model_highDelta)
model_sol_highDelta <- model_solve(model_highDelta)

Price.ZC <- varphi(model_sol_highDelta,
                   omega.varphi=omega_ZCB,
                   H = H/model_sol$tstep)
yds_CR_EZ_highDelta <- Price.ZC$r.t/100

EV   <- EV.fct(model_sol_highDelta,h=H/model_sol$tstep)
ET.P_CR_EZ <- EV$EX$T_at[1:(H/model_sol$tstep)]
ET.Q_CR_EZ <- varphi.tilde(model_sol_highDelta,omega_T.at,H/model_sol$tstep)[[1]]/
  varphi(model_sol_highDelta,omega_ZCB,H/model_sol$tstep)[[3]]
TRP_CR_EZ_highDelta <- ET.Q_CR_EZ - ET.P_CR_EZ

# CR model CRRA ----------------------------------------------------------------

gamma <- 1.01
gamma <- 1.45

model.CRRA <- model_sol
model.CRRA$parameters$gamma <- gamma
model.CRRA$parameters$delta <- model_sol$parameters$delta
model.CRRA <- solveParam4c(model.CRRA,indic_CRRA=TRUE)
model_CRRA_sol <- model_solve(model.CRRA,
                              indic_mitig = T,
                              indic_CRRA = T)
Price.ZC <- varphi(model_CRRA_sol,
                   omega.varphi=omega_ZCB,
                   H = H/model_sol$tstep)

yds_CR_CRRA <- Price.ZC$r.t/100

EV   <- EV.fct(model_CRRA_sol,h=H/model_sol$tstep)
ET.P_CR_CRRA <- EV$EX$T_at[1:(H/model_sol$tstep)]
ET.Q_CR_CRRA <- varphi.tilde(model_CRRA_sol,omega_T.at,H/model_sol$tstep)[[1]]/
  varphi(model_CRRA_sol,omega_ZCB,H/model_sol$tstep)[[3]]
TRP_CR_CRRA <- ET.Q_CR_CRRA - ET.P_CR_CRRA



# ------------------------------------------------------------------------------
# Plots ----
# 1----
FILE = "/outputs/Figures/Figure_YC_RF.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=11, width=11, height=5)

par(mfrow=c(1,2))

plot(1:H,rep(0,H),ylim=c(0,5),col="white",
     ylab="Interest rate (in percent)",
     xlab="Maturity, in years",las=1,
     main="(a) Term structure of real rates",cex.main=1.2)

polygon(c(seq(model_sol$tstep,H,by=model_sol$tstep),
          seq(H,model_sol$tstep,by=-model_sol$tstep)),
        100*c(yds_CR_EZ_lowDelta,rev(yds_CR_EZ_highDelta)),
        col="light grey",border = NaN)

# Bauer and Rudebusch yield curves:
BR_sdrs <- read.csv("data/Bauer_Rudebusch_sdrs.csv")
lines(BR_sdrs$uc_y10_2019[1:H],col=col.BR,lwd=2,lty=3)
lines(BR_sdrs$arbreak_y10_2019[1:H],col=col.BR,lwd=2,lty=2)
lines(BR_sdrs$arlearn_y10_2019[1:H],col=col.BR,lwd=2,lty=1)

lines(100*yds_lemoine,type="l",col=col.Lemoine,lwd=2)
lines(100*yds_BKO,col=col.BKO,lwd=2)
lines(100*yds_BKO_highDamages,col=col.BKO,lwd=2,lty=2)
lines(Price.ZC$date-model_sol$vec_date[1],100*yds_CR_EZ,col=col.CR,lwd=2)
lines(Price.ZC$date-model_sol$vec_date[1],100*yds_CR_CRRA,col=col.CR,lwd=2,lty=2)

legend("topright",
       legend=c(expression(paste("CR, Epstein-Zin (",gamma," = 7)",sep="")),
                expression(paste("CR, CRRA (",gamma," = 1.45)",sep="")),
                "Lemoine (2021)",
                "BKO (2019)",
                "BKO (2019), High dmg",
                "",
                "BR (2022), AR-learning",
                "BR (2022), AR-break",
                "BR (2022), UC"),
       lty=c(1:2,
             1,
             1:2,NaN,
             1:3),
       cex=1,
       col=c(col.CR,col.CR,
             col.Lemoine,
             col.BKO,col.BKO,"white",
             col.BR,col.BR,col.BR),
       lwd=2,
       bty = "n", bg="white",
       ncol=2)

grid()


# Make plots with expected path of of future 10-year yield:
EV <- EV.fct(model_sol)
Price.ZC2 <- varphi(model_sol,
                   omega.varphi=matrix(0,model_sol$n.X,1),
                   H = model_sol$horiz.2100)
maturities = c(2,10)
expected.yds <- compute_CMT(model_sol,EV,Price.ZC2,maturities)
make_figure_CMT(expected.yds,model_sol,maturities,
                indic_only_first = TRUE)

dev.off()


# 2----
FILE = "/outputs/Figures/Figure_TRP_comparison.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=11, width=7, height=4)

par(mfrow=c(1,1))
par(plt=c(.15,.95,.2,.95))

# Temperature risk premiums:

ylim <- c(0,max(TRP_BKO,TRP_CR_CRRA,TRP_CR_EZ,TRP_lemoine))
plot(1:H,rep(0,H),ylim=ylim,col="white",
     ylab="Temperature risk premium, in °C",
     xlab="Maturity, in years",las=1)

polygon(c(seq(model_sol$tstep,H,by=model_sol$tstep),
          seq(H,model_sol$tstep,by=-model_sol$tstep)),
        c(TRP_CR_EZ_lowDelta,rev(TRP_CR_EZ_highDelta)),
        col="light grey",border = NaN)


lines(TRP_lemoine,type="l",col=col.Lemoine,lwd=2)
lines(TRP_BKO,col=col.BKO,lwd=2)
lines(TRP_BKO_highDamages,col=col.BKO,lwd=2,lty=2)

lines(Price.ZC$date-model_sol$vec_date[1],TRP_CR_EZ,
      col=col.CR,lwd=2)
lines(Price.ZC$date-model_sol$vec_date[1],TRP_CR_CRRA,
      col=col.CR,lwd=2,lty=2)

grid()


legend("topleft",
       legend=c(expression(paste("CR, Epstein-Zin (",gamma," = 7)",sep="")),
                expression(paste("CR, CRRA (",gamma," = 1.45)",sep="")),
                "Lemoine (2021)",
                "BKO (2019)",
                "BKO (2019), High damages"),
       lty=c(1:2,
             1,
             1:2),
       cex=1,
       col=c(col.CR,col.CR,
             col.Lemoine,
             col.BKO,col.BKO),
       lwd=2,
       bty = "n", bg="white",
       ncol=1)

dev.off()
