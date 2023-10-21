
#model$parameters$gamma   <- 10
#model$parameters$mu_d    <- .02
#model$parameters$sigma_a <- .05


# # Check influence of mu_D on long-term rates:
# all.mu <- seq(0,10,by=,1)
# res.pricing <- varphi(model_sol,omega_ZCB,H=20)
# plot(res.pricing$r.t,ylim=c(-20,5))
# for(mu_n in all.mu){
#   print(mu_n)
#   model_new <- model
#   model_new$parameters$mu_n <- mu_n
#   model_new_sol<-model_solve(model_new,theta0)
#   
#   res.pricing <- varphi(model_new_sol,omega_ZCB,H=20)
#   lines(res.pricing$r.t)
#   
#   print(scc.fct(model_new_sol,h=0))
# }

# # Check that E^Q(T)=E^P(T) when no risk premiums:
# H <- 100
# model_test <- model
# model_test$parameters$mu_d <- 0
# model_test$parameters$mu_n <- 0
# model_test$parameters$q0   <- 0
# model_test$parameters$b_sk <- 0
# model_test_sol<-model_solve(model_test,theta0)
# ETQ <- varphi.tilde(model_test_sol,omega_T.at,H=H)[[1]]/varphi(model_test_sol,0*omega_T.at,H=H)[[3]]
# EV  <- EV.fct(model_test_sol,H)
# ETP <- EV$EX$T_at
# plot(ETP)
# lines(ETQ)


# Inflation linked bonds:

#model$parameters$mu_d <- .015

mu_PI <- list(muprice_0 = 0.015*model$tstep,
              muprice_1 = matrix(0,model_sol$n.X,1))
indic_D <- which(model$names.var.X=="D")
mu_PI$muprice_1[indic_D] <- 0.1
mu_PI$muprice_1[model_sol$n.Z+1] <- -.01*0 # eta_A
model_sol<-model_solve(model,theta0,mu_altern = mu_PI)

omega_PI <- matrix(0,model_sol$n.X,1)
omega_PI[which(model$names.var.X=="Cum_dc")] <- 1 # where PI_t is.
H <- 20
prices.nom.bonds <- varphi(model_sol,omega.varphi = - omega_PI,H=H)
nom.rates <- prices.nom.bonds$r.t - 100/((1:H)*model$tstep)*c(t(mu_PI$muprice_1)%*%model_sol$X)

prices.ZCRF.bonds <- varphi(model_sol,omega.varphi = omega_ZCB,H=H)

plot(nom.rates,ylim=c(min(prices.ZCRF.bonds$r.t),
                      max(nom.rates)))
lines(prices.ZCRF.bonds$r.t)

plot(nom.rates - prices.ZCRF.bonds$r.t)

#Model-implied EV
EV<-EV.fct(model_sol,H)
# Get matrix of forecasts:
EXh <- matrix(unlist(EV$EXh),model_sol$n.X,H)
expected.yearly.inflation <- 100*(mu_PI$muprice_0 +
                                    matrix(mu_PI$muprice_1,nrow = 1) %*% EXh)/model$tstep
lines(c(expected.yearly.inflation))


stop()





# Attempt of Merton model:

# First step: compute the sequence of mu_A0t's:
mu_A <- list(muprice_0 = 0,
             muprice_1 = matrix(0,model_sol$n.X,1))
indic_D <- model_sol$n.Z + model_sol$n.W - 1

mu_A$muprice_1[1] <- 3
#mu_A$muprice_1[indic_D] <- - 4

model_sol<-model_solve(model,theta0,mu_altern = mu_A)
#EV<-EV.fct(model_sol,H)

omega_A <- matrix(0,model_sol$n.X,1)
omega_A[13] <- 1 # where A_t is (tilde{A} at that stage).
H <- 20
prices.Atilde <- varphi(model_sol,omega.varphi = omega_A,H=H)
varphi.tilde.A   <-  prices.Atilde$P.t
varphi.tilde.A_1 <- c(1,varphi.tilde.A[1:(H-1)])
mu_A0t <- - log(varphi.tilde.A) + log(varphi.tilde.A_1)
# Update mu_A0t:
mu_A$muprice_0 <- c(0,mu_A0t,rep(mu_A0t[H],model$Tmax-H-1))
model_sol<-model_solve(model,theta0,mu_altern = mu_A)
EV<-EV.fct(model_sol,H)

# Check Euler equations:
prices.Atilde <- varphi(model_sol,omega.varphi = omega_A,H=H)

# Compare expected trajectories of exp(hr_{t,h}) and A:
expected.A <- NULL
for(h in 1:H){
  expected.A.h <- multi.lt.fct.N(model_sol,U=omega_A,h=h)
  expected.A <- c(expected.A,expected.A.h)
}
plot(1/prices.ZCRF.bonds$P.t,ylim=c(min(1/prices.ZCRF.bonds$P.t,expected.A),
                                    max(1/prices.ZCRF.bonds$P.t,expected.A)))
lines(expected.A)

# Price defaultable bonds:
A_bar <- seq(.5,3,length.out=10)
x <- exp(seq(-20,30,length.out = 4000))                                         
x <- exp(seq(-10,20,length.out = 10000))                                         
x <- exp(seq(-7,12,length.out = 3000))
H <- 10
yy <- varphi(model_sol,omega_ZCB,H=H)
B.star <- matrix(NaN,H,length(A_bar))
for(h in 1:H){
  print(h)
  xx.fast <- varphi.hat.fast(model_sol,omega = omega_ZCB,
                             H = h,x = x,
                             a = - omega_A,b = - log(A_bar))
  xx2.fast <- varphi.hat.fast(model_sol,omega = omega_A,
                              H = h,x = x,
                              a = omega_A,b = log(A_bar))
  
  B.star[h,] <- xx.fast + 1/A_bar * xx2.fast
}
plot(c(NaN,yy$P.t[2:H]),ylim=c(min(B.star,yy$P.t[2:H],na.rm = TRUE),
                       max(B.star,yy$P.t[2:H],na.rm = TRUE)))
lines(B.star[,1])
for(i in 1:length(A_bar)){
  lines(B.star[,i],col=i)
}


# Price equity:
xx.fast.E <- varphi.hat.fast(model_sol,omega = omega_A,
                             H = h,x = x,
                             a = - omega_A,b = - log(A_bar))
E <- xx.fast.E - A_bar*xx.fast

# Expected value of equity at maturity:
# fourier(model_sol,) 
EA1 <- fourier_complete(model_sol,x,omega_A,-omega_A,-log(A_bar),h=H)
EA2 <- fourier_complete(model_sol,x,0*omega_A,-omega_A,-log(A_bar),h=H)
ExpEqH <- EA1 - A_bar*EA2
# Expected return on Equity:
return.Eq.yearly <- log(ExpEqH/E)/(model$tstep*H)
plot(A_bar,return.Eq.yearly)



