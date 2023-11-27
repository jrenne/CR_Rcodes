

model_sol_new <- model_sol
model_sol_new$parameters$b_sk <- .000001

# model_sol_new$target_vector["ECumD2"] <- .995
# model_sol_new$target_vector["ECumD4"] <- .99
# model_sol_new$target_vector["stdCumD4"] <- .001

model_sol_new <- solveParam4D(model_sol_new)
model_sol_new <- model_solve(model_sol_new)



H <- 100

# Returns of risk-free strategy:
omega.ZC <- matrix(0,model_sol_new$n.X,1)
prices.ZCRF.bonds <- varphi(model_sol_new,
                            omega.varphi = omega.ZC,
                            H = H)

plot(prices.ZCRF.bonds$P.t,type="l")


# Computation of wealth:

# Selection vector for Cum_dc:
omega_C <- matrix(0,model_sol_new$n.X,1)
omega_C[which(model_sol_new$names.var.X=="Cum_dc")] <- 1
#omega_C[which(model_sol_new$names.var.X=="T_at")] <- 1

X.shock <- model_sol_new$X
X.shock[which(model$names.var.X=="E")] <- - 1 + X.shock[which(model$names.var.X=="E")]

prices.C   <- varphi(model_sol_new,
                     omega.varphi = omega_C,
                     H = H)
prices.C.shock   <- varphi(model_sol_new,
                           omega.varphi = omega_C,
                           H = H,
                           X = X.shock)
D <- sum(prices.C$P.t) - sum(prices.C.shock$P.t)

gamma <- model_sol$parameters$gamma
delta <- model_sol$parameters$delta
EC       <- NULL
EC.shock <- NULL
for(h in 1:H){
  Uh <- matrix(model_sol$mu_c1,model_sol$n.X,h)
  res.lt       <- multi.lt.fct.Uh(model_sol,Uh,X=model_sol_new$X,t=0)
  res.lt.shock <- multi.lt.fct.Uh(model_sol,Uh,X=X.shock,t=0)
  EC <- c(EC,res.lt$uX_t.h)
  EC.shock <- c(EC.shock,res.lt.shock$uX_t.h)
}
#cbind(EC,EC.shock)

stop()




nb.traj  <- 10000
nb.simul <- H
test <- simul.function(model_sol_new,
                       nb.simul,
                       nb.traj)
model_sol_new_shock <- model_sol_new
model_sol_new_shock$X <- X.shock
test.shock <- simul.function(model_sol_new_shock,
                             nb.simul,
                             nb.traj)

plot(density(test$Cum_dc[H,]))
lines(density(test.shock$Cum_dc[H,]),col="red")

print(cbind(apply(test$Cum_dc,1,mean),
            apply(test.shock$Cum_dc,1,mean)))
plot(apply(test$Cum_dc,1,mean),type="l")
lines(apply(test.shock$Cum_dc,1,mean),col="red")

EC.shock <- apply(test.shock$Cum_dc,1,mean)
EC       <- apply(test$Cum_dc,1,mean)
ExpBenef <- EC.shock - EC
plot(ExpBenef,type="l")

DF <- varphi(model_sol_new,
             omega.varphi = omega.ZC,
             H = H)$P.t
DF.shock <- varphi(model_sol_new_shock,
                   omega.varphi = omega.ZC,
                   H = H)$P.t

sum(DF.shock*EC.shock - DF*EC)




