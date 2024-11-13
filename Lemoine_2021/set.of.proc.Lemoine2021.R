
simul.model.Lemoine <- function(model,X0,H,nb.replic=1,damage.type="Lemoine"){
  
  all.C     <- matrix(NaN,H,nb.replic)
  all.Cexog <- matrix(NaN,H,nb.replic)
  all.L  <- matrix(NaN,H,nb.replic)
  all.T  <- matrix(NaN,H,nb.replic)
  all.M1 <- matrix(NaN,H,nb.replic)
  all.M2 <- matrix(NaN,H,nb.replic)
  
  C     <- matrix(X0[1],1,nb.replic)
  Cexog <- matrix(X0[1],1,nb.replic)
  L  <- matrix(X0[2],1,nb.replic)
  T  <- matrix(X0[3],1,nb.replic)
  M1 <- matrix(X0[4],1,nb.replic)
  M2 <- matrix(X0[5],1,nb.replic)
  
  mut      <- model$mu0    * exp(- model$mu1    * (1:H + 2014 - 1960))
  gammat   <- model$gamma0 * exp(- model$gamma1 * (1:H + 2014 - 1960))
  delta.lt <- model$l0     * exp(- model$l1     * (1:H + 2014 - 1960))
  
  for(t in 1:H){
    if(damage.type=="Lemoine"){
      C <- C * (1 + mut[t] - model$alpha * (T - T0) + model$sigma * rnorm(nb.replic))
      Cexog <- C
    }
    if(damage.type=="G-N"){
      C <- C * (1 + mut[t] - .00026 * T + model$sigma * rnorm(nb.replic))
      Cexog <- C
    }
    if(damage.type=="G-DJO"){
      C <- C * (1 + mut[t] - .00137 * T + model$sigma * rnorm(nb.replic))
      Cexog <- C
    }
    if(damage.type=="G-W"){
      C <- C * (1 + mut[t] - .000075 * T^3.25 + model$sigma * rnorm(nb.replic))
      Cexog <- C
    }
    if(damage.type=="L-N"){
      Cexog <- Cexog * (1 + mut[t] + model$sigma * rnorm(nb.replic))
      C <- Cexog * 1/(1+.00266*T^2)
    }
    if(damage.type=="L-W"){
      Cexog <- Cexog * (1 + mut[t] + model$sigma * rnorm(nb.replic))
      C <- Cexog * 1/(1+(T/20.64)^2+(T/6.081)^6.754)
    }
    if(damage.type=="L-Barnett"){
      Cexog <- Cexog * (1 + mut[t] + model$sigma * rnorm(nb.replic))
      y.bar <- 2
      C <- Cexog * exp(-0.0127*T -0.0005*T^2 - .02 * (T - y.bar)^2 * (T > y.bar))
    }
    if(damage.type=="L-Bansal"){
      Cexog <- Cexog * (1 + mut[t] + model$sigma * rnorm(nb.replic))
      pi.t <- .0050 + .0033 * T
      N <- rpois(nb.replic,pi.t)
      d.t  <- (.0011 * T + .0011 * T^2) * (T > 2)
      D <- rgamma(nb.replic,shape=1,scale=d.t*N)
      C <- Cexog * exp(-D)
    }
    L <- L * exp(delta.lt[t])
    M1 <- M1 + model$psi1 * gammat[t] * C
    M2 <- M2 + model$psi2 * (1 - model$psi1) * gammat[t] * C -
      model$psi * M2
    M  <- M1 + M2
    #F  <- model$nu * (log(model$M0/model$Mpre) + M/model$M0 - 1)
    F  <- model$nu * log(M/model$Mpre)
    T  <- T + model$phi * (F - model$s * T)
    
    all.C[t,]      <- C
    all.Cexog[t,]  <- Cexog
    all.L[t,]  <- L
    all.T[t,]  <- T
    all.M1[t,] <- M1
    all.M2[t,] <- M2
  }
  return(list(
    all.C  = all.C,
    all.L  = all.L,
    all.T  = all.T,
    all.M1 = all.M1,
    all.M2 = all.M2
  ))
}

make.model.Lemoine <- function(indic.uncertainty.alpha = 1,
                               indic.uncertainty.S = 1,
                               indic.uncertainty.C = 1,
                               nb_discret_values = 6,
                               H = 200,
                               nb.replic = 500, setseed = NaN){
  
  damage.type <- "Lemoine"
  # damage.type <- "G-N"
  # damage.type <- "G-W"
  # damage.type <- "L-N"
  # damage.type <- "L-W"
  # damage.type <- "G-DJO"
  # damage.type <- "L-Barnett"
  # damage.type <- "L-Bansal"
  
  mu0    <- .047
  mu1    <- .0134
  alpha  <- .0018
  sigma  <- .014 * indic.uncertainty.C
  l0     <- .022
  l1     <- .012
  psi1   <- .2
  psi2   <- .393
  psi    <- .0023
  nu     <- 5.35
  Mpre   <- 605
  gamma0 <- .28
  gamma1 <- .015
  phi    <- .0394
  S      <- 3.13
  
  #M0     <- 800 # where radiative forcings are linearized
  
  s <- nu * log(2) / S
  
  C0  <- 55 * 1.1772 # in trillion USD 2014
  L0  <- 7.1 * 10^9
  T0  <- .8794
  M10 <- 706
  M20 <- 148
  # NB: M = M10 + M20, and divide by 2.13 to have M expressed in ppmv CO2
  
  eta <- 1.45
  beta <- 1 - .015
  
  X0 <- c(C0,L0,T0,M10,M20)
  
  model <- list(
    mu0 = mu0,
    mu1 = mu1,
    alpha = alpha,
    sigma = sigma,
    l0 = l0,
    l1 = l1,
    psi1 = psi1,
    psi2 = psi2,
    psi = psi,
    nu = nu,
    Mpre = Mpre,
    gamma0 = gamma0,
    gamma1 = gamma1,
    phi = phi,
    #M0 = M0,
    s = s
  )
  
  if(!is.na(setseed)){
    set.seed(setseed)
  }
  res.simul <- simul.model.Lemoine(model,X0,H=H,nb.replic = nb.replic,
                                   damage.type = damage.type)
  
  res.mean.traj <- lapply(res.simul,function(x){apply(x,1,mean)})
  
  return(list(res.simul=res.simul,
              res.mean.traj,res.mean.traj))
}

