# tempo check conditional means:

# model$parameters$mu_T <- .00000000000001
# model$parameters$mu_D <- .00000000000001
# model$parameters$mu_H <- .00000000000001
# model$parameters$mu_N <- .00000000000001
# model$parameters$a_N <- .0000001
# model$parameters$b_N <- 0*model$parameters$b_N
# model$parameters$sigma_a <- .00000000000001
model_sol <- model_solve(model,indic_CRRA = FALSE)

nb.t <- 100
nb.traj <- 1000

res <- simul.function(model_sol,nb.simul.t=nb.t,nb.traj=nb.traj)#,setseed=123)

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

EV <- EV.fct(model_sol,h=nb.t)
lines(EV$EX$T_at)

EV_shock <- EV.fct(model_sol_shock,h=nb.t)
lines(ETAT_shock.sim,col="red")
lines(EV_shock$EX$T_at,col="red")

plot(ETAT_shock.sim - ETAT.sim)
lines(EV_shock$EX$T_at-EV$EX$T_at)

