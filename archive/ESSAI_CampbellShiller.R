

N <- 100

tau_bar <- log(8)

sigma <- .3

tau <- tau_bar + sigma * rnorm(N)

kappa1 <- exp(tau_bar)/(1 + exp(tau_bar))
kappa0 <- log(1+exp(tau_bar)) - kappa1 * tau_bar

P.over.D <- exp(tau)
plot(P.over.D * 5)

P.over.D_tp1 <- P.over.D[2:N]
P.over.D_t   <- P.over.D[1:(N-1)]

r <- (1 + P.over.D_tp1)/P.over.D_t - 1
r_approx <- kappa0 + kappa1*log(P.over.D_tp1) - log(P.over.D_t)

plot(r,type="l",ylim=c(min(r,r_approx),max(r,r_approx)))
lines(r_approx,col="red")
