
gamma <- 1.24
rho <- .002

# Nordhaus 2017b: D = .00236*T^2

# gamma <- 0.001
# rho <- .03

# gamma <- 1.57
# rho <- .008

#gamma <- 7
model.CRRA <- model
model.CRRA$parameters$gamma <- gamma
model.CRRA$parameters$delta      = (1 - rho)^5

model.CRRA$target_vector["mu_c0"] <- .015*model$tstep
#model.CRRA$target_vector["sigma_c0"] <- sqrt(3*.014)
#model.CRRA$target_vector["sigma_c0"] <- sqrt(5*.014)
model.CRRA <- solveParam4c(model.CRRA,indic_CRRA=TRUE)

# model.CRRA$parameters$a_D  <- gamma * model.CRRA$parameters$a_D * 0
# model.CRRA$parameters$b_D  <- gamma * .0018 * model.CRRA$tstep
# model.CRRA$parameters$mu_D <- gamma * model.CRRA$parameters$mu_D * 4
# #model$parameters$mu_T <- model$parameters$mu_T * 3
model.CRRA$parameters$a_D <- model$parameters$a_D/2
model.CRRA$parameters$b_D <- model$parameters$b_D/2
model.CRRA$parameters$mu_D <- .0001
model.CRRA$parameters$b_sk <- 0
# model.CRRA$parameters$a_N <- 0.0001
# model.CRRA$parameters$b_N <- 0.0001

H <- (2300 - 2020)/model$tstep


# Determine mitigation rate to have appropriate dynamics of Emissions & T_AT:
theta0  <- c(log(0.17),-1/21*log(0.17))
#theta0  <- c(log(0.17),-1/10*log(0.17))
#theta0  <- c(log(0.17),-1/5*log(0.17))
a <- abs(theta0[1])
b <- abs(theta0[2])
mu.chosen <- pmin(exp(- a + b*(1:model_sol$Tmax)),1)
plot(mu.chosen)

# Solve model:
model_CRRA_sol <- model_solve(model.CRRA,
                              indic_mitig = FALSE,
                              indic_CRRA = TRUE,
                              mu.chosen = mu.chosen)

EV <- EV.fct(model_sol = model_CRRA_sol,h = H)
par(mfrow=c(2,3))
plot(EV$EX$E_ind)
plot(EV$EX$T_at)
plot(cumsum(EV$EX$D))

res.scc <- scc.fct.CRRA(model_CRRA_sol,t = 0,H = H)
print(res.scc$SCC.CO2)
plot(res.scc$scc.decomp)
omega.ZC <- matrix(0,model_CRRA_sol$n.X,1)
prices.ZCRF.bonds   <- varphi(model_CRRA_sol,
                              omega.varphi = omega.ZC,
                              H = H)
plot(prices.ZCRF.bonds$r.t)
plot(1:20,prices.ZCRF.bonds$r.t[1:20])


