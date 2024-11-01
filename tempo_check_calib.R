


T0    <- .5
Tstar <- 3
tstar <- 16
kappa <- .9

alpha <- .2

Tinf <- (Tstar - T0 * exp(-alpha*tstar))/(1 - exp(-alpha*tstar))

seq_Ti_0_tstar <- Tinf - (Tinf - T0)*exp(-alpha*0:tstar)


# For D:

tstar*Tinf - (Tinf - T0)*(1 - exp(-alpha*tstar))/(1 - exp(-alpha))

linearseqT <- seq(T0,Tstar,length.out=tstar+1)
sum(linearseqT[1:tstar])

(Tstar*(tstar-1) + T0*(tstar+1))/2


a <- (tstar*(1 - exp(-alpha)) - (1- exp(-alpha*tstar)))/
  ((1- exp(-alpha*tstar))*(1- exp(-alpha)))

b <- ((1 - exp(-alpha*tstar)) - tstar*exp(-alpha*tstar)*(1-exp(-alpha)))/
  ((1- exp(-alpha*tstar))*(1- exp(-alpha)))

a*Tstar + b*T0
sum(seq_Ti_0_tstar[1:tstar])


# For N:

t <- tstar + 10
seq_Ti_0_t <- Tinf - (Tinf - T0)*exp(-alpha*0:t)

a <- 1/(1 - exp(-alpha*tstar))*(
  (1 - kappa^t)/(1 - kappa) -
    (1 - kappa^t*exp(-alpha*t))/(1 - kappa*exp(-alpha))
)
b <- - 1/(1 - exp(-alpha*tstar))*(
  (1 - kappa^t)*exp(-alpha*tstar)/(1 - kappa) -
    (1 - kappa^t*exp(-alpha*t))/(1 - kappa*exp(-alpha))
)

a*Tstar + b*T0
sum(seq_Ti_0_t[1:t]*kappa^(0:(t-1)))



