

model_sol_new <- model_sol
#model_sol_new$parameters$b_sk <- .000001

model_sol_new$parameters$mu_T <- .0000001
model_sol_new$parameters$mu_N <- .0000001
model_sol_new$parameters$mu_H <- .0000001
model_sol_new$parameters$mu_D <- .05
model_sol_new$parameters$sigma_a <- .0000001

gamma <- 2
model_sol_new$parameters$gamma          <- gamma
model_sol_new$target_vector["mu_c0"]    <- .06
model_sol_new$target_vector["sigma_c0"] <- .00001
model_sol_new <- solveParam4c(model_sol_new,
                              indic_CRRA=TRUE)
model_sol_new$parameters$a_D  <- gamma * model_sol_new$parameters$a_D * 0
model_sol_new$parameters$b_D  <- gamma * .0018 * model$tstep
# model_sol_new$parameters$a_N  <- 0
# model_sol_new$parameters$b_N  <- 0
model_sol_new$parameters$b_sk <- 0
model_sol_new <- model_solve(model_sol_new,
                             indic_mitig = T,
                             indic_CRRA = T)

#model_sol_new <- model_solve(model_sol_new)

H <- 100

EV<-EV.fct(model_sol_new,H)
par(mfrow=c(2,3))
plot(model$vec_date[1:H],EV$EX$T_at,type="l")
lines(model$vec_date[1:H],EV$EX$T_atW,type="l",col="red")
lines(model$vec_date[1:H],EV$EX$T_at+2*sqrt(EV$VX$T_at),type="l",col="red",lty=2)
lines(model$vec_date[1:H],EV$EX$T_at-2*sqrt(EV$VX$T_at),type="l",col="red",lty=2)
plot(model$vec_date[1:H],EV$EX$delc,type="l",ylim=c(-.02,.07))
lines(model$vec_date[1:H],EV$EX$delc+2*sqrt(EV$VX$delc),type="l",col="red",lty=2)
lines(model$vec_date[1:H],EV$EX$delc-2*sqrt(EV$VX$delc),type="l",col="red",lty=2)
plot(model$vec_date[1:H],EV$EX$Cum_D,type="l")
lines(model$vec_date[1:H],EV$EX$Cum_D+2*sqrt(EV$VX$Cum_D),type="l",col="red",lty=2)
lines(model$vec_date[1:H],EV$EX$Cum_D-2*sqrt(EV$VX$Cum_D),type="l",col="red",lty=2)
plot(model$vec_date[1:H],EV$EX$H,type="l")
lines(model$vec_date[1:H],EV$EX$H+2*sqrt(EV$VX$H),type="l",col="red",lty=2)
lines(model$vec_date[1:H],EV$EX$H-2*sqrt(EV$VX$H),type="l",col="red",lty=2)


#stop()



# Compute prices of risk-free bonds:
omega.ZC <- matrix(0,model_sol_new$n.X,1)
prices.ZCRF.bonds <- varphi(model_sol_new,
                            omega.varphi = omega.ZC,
                            H = H)
plot(prices.ZCRF.bonds$P.t,type="l")

# Selection vector for Cum_dc:
omega_C <- matrix(0,model_sol_new$n.X,1)
omega_C[which(model_sol_new$names.var.X=="Cum_dc")] <- 1
#omega_C[which(model_sol_new$names.var.X=="T_at")] <- 1

epsilon <- .01
X.shock <- model_sol_new$X
X.shock[which(model$names.var.X=="M_at")] <- - epsilon + 
  X.shock[which(model$names.var.X=="M_at")]
prices.C   <- varphi(model_sol_new,
                     omega.varphi = omega_C,
                     H = H)
prices.C.shock   <- varphi(model_sol_new,
                           omega.varphi = omega_C,
                           H = H,
                           X = X.shock)
D <- sum(prices.C$P.t) - sum(prices.C.shock$P.t)

EC       <- NULL
EC.shock <- NULL
for(h in 1:H){
  Uh <- matrix(model_sol_new$mu_c1,model_sol$n.X,h)
  res.lt       <- multi.lt.fct.Uh(model_sol_new,Uh,X=model_sol_new$X,t=0)
  res.lt.shock <- multi.lt.fct.Uh(model_sol_new,Uh,X=X.shock,t=0)
  EC <- c(EC,res.lt$uX_t.h)
  EC.shock <- c(EC.shock,res.lt.shock$uX_t.h)
}
#cbind(EC,EC.shock)


# Discounted expected benefits:
NPV <- sum(prices.ZCRF.bonds$P.t * (EC.shock - EC)/epsilon) *
  10^3 * model$parameters$c0 / 3.667
SCC <- scc.fct.CRRA(model_sol_new,H=H)
cbind(SCC$SCC.CO2,NPV)

plot(SCC$scc.decomp / 3.667)
lines(prices.ZCRF.bonds$P.t * (EC.shock - EC)/epsilon *
        10^3 * model$parameters$c0 / 3.667)

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




