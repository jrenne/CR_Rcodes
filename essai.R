# ==============================================================================
# Check sensitivity to assumptions regarding mu_t:
# (i) parametric function + (ii) time consistency
# ==============================================================================


EV <- EV.fct(model_sol)

t <- 10 # t=0 corresponds to initial date


# Note: the first entry of mu is for t=0 (it is replaced with mu0)

if(t>0){
  Xt <- EV$EXh[[t]]
}else{
  Xt <- model_sol$X
}
Xt[6] <- 10 * Xt[6]
#Xt[10] <- 10*Xt[10]

RES.2param <- res.optim(model_sol,
                        model_sol$theta0,
                        Tend = model_sol$Tmax - t,
                        X = Xt)

# use of a logit function to code mu_t:
theta <- .5 + .99*(model_sol$mu[(t+1):model_sol$Tmax] - .5)
theta <- log(theta/(1-theta))
RES.Tmax.param <- res.optim(model_sol,
                            theta,
                            Tend = model_sol$Tmax - t,
                            X = Xt)

# theta <- RES.Tmax.param$par
# theta[1] <- 10
# utility.optim(model_sol,
#               theta,
#               Tend = model_sol$Tmax - t,
#               X    = Xt)

# mu_aux <- mu.function(model_sol,
#                       theta = theta,
#                       t.ini = t)
# essai <- mu_dep(model_sol,mu_aux)
# plot(extract(essai$omega0,1))

plot(mu.function(model_sol,
                 theta = RES.2param$par,
                 t.ini = t))
lines(mu.function(model_sol,
                  theta = RES.Tmax.param$par,
                  t.ini = t))


stop()

cbind(RES$par,model_sol$theta.opt)

print(RES$value)

plot(mu.function(model_sol,
                 theta = model_sol$theta.opt,
                 t.ini = t))
lines(mu.function(model_sol,
                  theta = RES$par,
                  t.ini = t))

stop()


utility.optim(model_sol,
              theta,
              Tend = model_sol$Tmax - t,
              X    = Xt)
plot(mu.function(model_sol,
                 theta,
                 t.ini = t + 5))
