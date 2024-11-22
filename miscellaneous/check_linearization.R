# ==============================================================================
# This script checks the quality of the approximation concerning the
# link between emissions and production.
# ==============================================================================

nb_replic <- 10

param <- model_sol$parameters
Tmax  <- model_sol$Tmax
tstep <- model_sol$tstep

gsigma <- matrix(0,Tmax,1)
sigma  <- matrix(0,Tmax,1)

#Carbon intensity
gsigma[1] <- param$gsigma1                                                      
for(i in 2:Tmax) gsigma[i] <- gsigma[i-1]*(1+param$delsigma)   
sigma[1] <- param$sigma0                                                         
for(i in 2:Tmax) sigma[i] <- sigma[i-1] * exp(gsigma[i-1] * tstep)      

sigma <- matrix(sigma[2:Tmax],ncol=1)

plot(sigma)

model_sol$mu_c

mu_ct    <- matrix(extract(model_sol$omega0,1),ncol=1)
sigma_ct <- matrix(extract(model_sol$omega,1),ncol=1)

vec_1 <- matrix(1,1,nb_replic)

E_ind <- matrix(NaN,Tmax,nb_replic)

# Exact simulation:
eps <- matrix(rnorm((Tmax-1)*nb_replic),Tmax-1,nb_replic)
mu  <- matrix(model_sol$mu[2:Tmax],ncol=1)
E_ind <- (sigma %*% vec_1) * (1 -  mu %*% vec_1) * 
  exp(apply((mu_ct %*% vec_1) + (sigma_ct %*% vec_1) * eps,2,cumsum))

# Approximate simulation:
lambda_t <- (sigma %*% vec_1) * (1 -  mu %*% vec_1) * 
  exp( apply((mu_ct %*% vec_1) + (sigma_ct^2 %*% vec_1)/2,2,cumsum))
y_tilde <- apply((sigma_ct %*% vec_1) * eps,2,cumsum)
E_ind_approx <- lambda_t * (1 + y_tilde)


par(mfrow=c(2,2))
for(i in 1:4){
  plot(E_ind[,i],type="l")
  points(E_ind_approx[,i],col="red")
}

        