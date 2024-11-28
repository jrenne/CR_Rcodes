# tempo check conditional means:

# model$parameters$mu_T <- .00000000000001
# model$parameters$mu_D <- .00000000000001
# model$parameters$mu_H <- .00000000000001
# model$parameters$mu_N <- .00000000000001
# model$parameters$a_N <- .0000001
# model$parameters$b_N <- 0*model$parameters$b_N
# model$parameters$sigma_a <- .00000000000001

#model_sol <- model_solve(model,indic_CRRA = FALSE)

nb.t <- 100
nb.traj <- 1000

res <- simul.function(model_sol,nb.simul.t=nb.t,nb.traj=nb.traj)#,setseed=123)
EV <- EV.fct(model_sol,h=nb.t)

plot(res$T_at[,10])
lines(EV$EX$T_at)

shock <- 1 # Gt Carbon

model_sol_shock <- model_sol
model_sol_shock$vector.ini$ini_Mat <- model_sol_shock$vector.ini$ini_Mat + shock
model_sol_shock$X[indic.M_at]      <- model_sol_shock$vector.ini$ini_Mat

model_sol_shock$vector.ini$ini_F <- model_sol_shock$vector.ini$ini_F - 
  model_sol$A0.star.inf[5,6] * shock
model_sol_shock$X[indic.Forc] <- model_sol_shock$vector.ini$ini_F

res_shock <- simul.function(model_sol_shock,nb.simul.t=nb.t,nb.traj=nb.traj)#,setseed=123)

diffTat <- res_shock$T_at - res$T_at

ETAT.sim <- apply(res$T_at,1,mean)
ETAT_shock.sim <- apply(res_shock$T_at,1,mean)

plot(ETAT.sim,type="l")

EV_shock <- EV.fct(model_sol_shock,h=nb.t)
lines(ETAT_shock.sim,col="red")
lines(EV_shock$EX$T_at,col="red")

plot(ETAT_shock.sim - ETAT.sim)
lines(EV_shock$EX$T_at-EV$EX$T_at)



# Check that SCC equal to sum of expected benefits when no uncertainty:

model_new <- model
model_new$parameters$mu_T <- .00000000000001
model_new$parameters$mu_D <- .00000000000001
model_new$parameters$mu_H <- .00000000000001
model_new$parameters$mu_N <- .00000000000001
model_new$parameters$a_N <- .0000001
model_new$parameters$b_N <- 0*model$parameters$b_N
model_new$parameters$sigma_a <- .00000000000001
model_sol_new <- model_solve(model_new,indic_CRRA = FALSE)

HH <- 300 # maximum horizon considered in infinite sums
epsilon <- 1
X.shock <- model_sol_new$X
X.shock[which(model_sol$names.var.X=="M_at")] <- - epsilon + 
  X.shock[which(model_sol$names.var.X=="M_at")]

omega_ZCB <- matrix(0,model$n.Z+model$n.W,1)
Omega_C   <- omega_ZCB
Omega_C[model_sol$names.var.X=="Cum_dc"] <- 1
prices.C   <- varphi(model_sol_new,
                     omega.varphi = Omega_C,
                     H = 300)
prices.C.shock   <- varphi(model_sol_new,
                           omega.varphi = Omega_C,
                           H = 300,
                           X = X.shock)
D <- sum(prices.C$P.t) - sum(prices.C.shock$P.t)

EC       <- NULL
EC.shock <- NULL
for(h in 1:HH){
  Uh <- matrix(model_sol_new$mu_c1,model_sol$n.X,h)
  res.lt       <- multi.lt.fct.Uh(model_sol_new,Uh,X=model_sol_new$X,t=0)
  res.lt.shock <- multi.lt.fct.Uh(model_sol_new,Uh,X=X.shock,t=0)
  EC <- c(EC,res.lt$uX_t.h)
  EC.shock <- c(EC.shock,res.lt.shock$uX_t.h)
}
# Discounted expected benefits:
omega_ZCB <- matrix(0,model$n.Z+model$n.W,1)
ZCB <- varphi(model_sol_new,omega_ZCB,HH)
NPV.CO2 <- sum(ZCB$P.t * (EC.shock - EC)/epsilon) *
  10^3 * model_sol$parameters$c0 / 3.667

