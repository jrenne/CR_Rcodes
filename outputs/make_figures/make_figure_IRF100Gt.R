# ==============================================================================
# Response of temperature to a 100-GtC shock (M_at)
# ==============================================================================

indic.M_at <- which(model_sol$names.var.X=="M_at")

shock <- 100

model_sol_shock <- model_sol
model_sol_shock$vector.ini$ini_Mat <- model_sol_shock$vector.ini$ini_Mat + shock
model_sol_shock$X[indic.M_at]      <- model_sol_shock$vector.ini$ini_Mat

EV       <- EV.fct(model_sol,h=50)
EV_shock <- EV.fct(model_sol_shock,h=50)

par(mfrow=c(2,2))
IRF <- EV_shock$EX$M_at - EV$EX$M_at
plot(IRF,type="l")
IRF <- EV_shock$EX$T_at - EV$EX$T_at
plot(IRF,type="l")


