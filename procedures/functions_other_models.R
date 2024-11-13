# ==============================================================================
# This script gathers functions associated with alternative modeling 
# frameworks, namely:
#
# Derek Lemoine, 2021.
# "The Climate Risk Premium: How Uncertainty Affects the Social Cost of Carbon,"
# Journal of the Association of Environmental and Resource Economists,
# University of Chicago Press, vol. 8(1), pages 27-57.
#
# Ravi Bansal & Marcelo Ochoa & Dana Kiku, 2019.
# "Climate Change Risk," 
# Work. Pap., Fuqua Sch. Bus., Duke Univ., Durham, NC.
# ==============================================================================


# BKO (2019) -------------------------------------------------------------------

make.BKO.model <- function(){
  # Baseline BKO model:
  
  model <- list(
    mu         = .015,
    sigma_eta  = .018,
    d          = -.05,
    ell0       = .01,
    ell1       = .01,
    nu         = .966,
    Theta      = 1,
    chi        = .2,
    sigma_zeta = 1,
    psi        = 1.5,
    gamma      = 5,
    delta      = .99
  )
  
  # Unconditional means:
  model$Ebar <- model$Theta * model$mu / (1 - model$nu)
  model$Tbar <- model$Ebar * model$chi
  
  return(model)
}

Psi.bko <- function(model,u,v){
  
  mu         <- model$mu
  sigma_eta  <- model$sigma_eta
  d          <- model$d
  ell0       <- model$ell0
  ell1       <- model$ell1
  nu         <- model$nu
  Theta      <- model$Theta
  chi        <- model$chi
  sigma_zeta <- model$sigma_zeta
  psi        <- model$psi
  gamma      <- model$gamma
  delta      <- model$delta
  
  a <- v*mu + u*chi*Theta*mu + (u*chi*Theta + v)^2/2*sigma_eta^2 +
    (u*chi)^2/2*sigma_zeta^2 + ell0*(exp(v*d) - 1)
  b <- u*nu + ell1*(exp(v*d) - 1)
  
  return(list(a=a,
              b=b))
}

phi.bko <- function(x){exp(x)-1}

solve_model.bko <- function(model){
  
  mu         <- model$mu
  sigma_eta  <- model$sigma_eta
  d          <- model$d
  ell0       <- model$ell0
  ell1       <- model$ell1
  nu         <- model$nu
  Theta      <- model$Theta
  chi        <- model$chi
  sigma_zeta <- model$sigma_zeta
  psi        <- model$psi
  gamma      <- model$gamma
  delta      <- model$delta
  
  Ebar <- model$Ebar
  Tbar <- model$Tbar
  
  theta <- (1 - gamma) / (1 - 1/psi)
  
  zbar <- 1
  for(i in 1:50){
    kappa1 <- exp(zbar)/(exp(zbar)-1)
    kappa0 <- kappa1 * zbar - log(exp(zbar)-1)
    
    A1 <- ell1/theta * phi.bko((1-gamma)*d)/(kappa1 - nu)
    A0 <- 1/(kappa1 - 1)*(
      log(delta) + kappa0 + (1 - 1/psi + chi*Theta*A1)*mu + ell0/theta*phi.bko((1-gamma)*d) +
        .5 * theta * ((1 - 1/psi + chi*Theta*A1)^2*sigma_eta^2 + (chi*A1*sigma_zeta)^2)
    )
    zbar <- A0 + A1 * Tbar
  }
  
  ab <- Psi.bko(model,(theta - 1)*A1, - gamma)
  a <- ab$a
  b <- ab$b
  
  rfbar <- - theta * log(delta) - (theta - 1)*kappa0 - 
    (theta - 1)*(1 - kappa1)*A0 - a +
    ((theta - 1)*kappa1*A1 - b) * Tbar
  
  a_m <- theta*log(delta) + (theta - 1)*(kappa0 + (1 - kappa1)*A0)
  b_m <- (theta - 1)*(-kappa1*A1)
  c_m <- (theta - 1)*A1
  d_m <- - gamma
  
  model_sol <- model
  
  model_sol$kappa0 = kappa0
  model_sol$kappa1 = kappa1
  model_sol$A0 = A0
  model_sol$A1 = A1
  model_sol$zbar = zbar
  model_sol$rfbar = rfbar
  model_sol$a_m = a_m
  model_sol$b_m = b_m
  model_sol$c_m = c_m
  model_sol$d_m = d_m
  
  return(model_sol)
}

compute_ab.bko <- function(model_sol,H,u=0){
  a <- 0
  b <- u
  all_a <- NULL
  all_b <- NULL
  
  a_m <- model_sol$a_m
  b_m <- model_sol$b_m
  c_m <- model_sol$c_m
  d_m <- model_sol$d_m
  
  for(h in 1:H){
    ab <- Psi.bko(model_sol,c_m+b,d_m)
    a <- a_m + a + ab$a
    b <- b_m + ab$b
    all_a <- c(all_a,a)
    all_b <- c(all_b,b)
  }
  return(list(
    all_a = all_a,
    all_b = all_b,
    all_r_a = - all_a/1:H,
    all_r_b = - all_b/1:H
  ))
}


simul.model.BKO <- function(model,
                            T0=model$Tbar,E0=model$Tbar,
                            H,nb.replic=1){
  
  all.delc <- matrix(NaN,H,nb.replic)
  all.C    <- matrix(NaN,H,nb.replic)
  all.D    <- matrix(NaN,H,nb.replic)
  all.T    <- matrix(NaN,H,nb.replic)
  all.E    <- matrix(NaN,H,nb.replic)
  
  T <- rep(T0,nb.replic)
  E <- rep(E0,nb.replic)
  C <- rep(1, nb.replic)
  
  for(t in 1:H){
    
    pit_1 <- model$ell0 + model$ell1*T
    D <- model$d * rpois(nb.replic,pmax(pit_1,0))
    eta  <- rnorm(nb.replic)
    zeta <- rnorm(nb.replic)
    
    delc <- model$mu + model$sigma_eta*eta + D
    C <- C*(1 + delc)
    
    E <- model$nu * E + model$Theta * (model$mu + eta) + 
      model$sigma_zeta * zeta
    T <- model$chi * E
    
    all.D[t,]    <- D
    all.delc[t,] <- delc
    all.C[t,]    <- C
    all.T[t,]    <- T
    all.E[t,]    <- E
  }
  return(list(all.delc = all.delc,
              all.C = all.C,
              all.D = all.D,
              all.T = all.T,
              all.E = all.E))
}

# Lemoine (2021)'s model -------------------------------------------------------

make.Lemoine.model <- function(){
  # Baseline model (Lemoine, 2021):
  
  model <- list(
    eta    = 1.45,
    beta   = 1 - .015,
    mu0    = .047,
    mu1    = .0134,
    sigma  = .014,
    l0     = .022,
    l1     = .012,
    psi1   = .2,
    psi2   = .393,
    psi    = .0023,
    nu     = 5.35,
    Mpre   = 605,
    gamma0 = .28,
    gamma1 = .015,
    phi    = .0394,
    S      = 3.13)
  
  model$s      <- model$nu * log(2) / model$S
  
  # Initial values:
  C0  <- 55 * 1.1772 # in trillion USD 2014
  L0  <- 7.1 * 10^9
  T0  <- .8794
  M10 <- 706
  M20 <- 148
  X0 <- c(C0,L0,T0,M10,M20)
  
  return(list(model=model,
              X0=X0))
}


simul.model.Lemoine <- function(model,
                                damage.type,
                                X0,H,nb.replic=1,
                                coef.multip = 1){
  
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
    
    D <- compute_alternative_damage(T,damage.type,
                                    time_step = 1,
                                    coef.multip = coef.multip)
    
    if(stringr::str_sub(damage.type,1,1)=="G"){
      # "Growth" specification
      C <- C * (1 + mut[t] - D + model$sigma * rnorm(nb.replic))
      Cexog <- C
    }else if(stringr::str_sub(damage.type,1,1)=="L"){
      # "Level" specification
      Cexog <- Cexog * (1 + mut[t] + model$sigma * rnorm(nb.replic))
      C <- Cexog * D
    }
    
    L <- L * exp(delta.lt[t])
    M1 <- M1 + model$psi1 * gammat[t] * C
    M2 <- M2 + model$psi2 * (1 - model$psi1) * gammat[t] * C -
      model$psi * M2
    M  <- M1 + M2
    # Linearized version:
    ##F  <- model$nu * (log(model$M0/model$Mpre) + M/model$M0 - 1)
    F  <- model$nu * log(M/model$Mpre)
    T  <- T + model$phi * (F - model$s * T)
    
    all.C[t,]     <- C
    all.Cexog[t,] <- Cexog
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

sumE <- function(x,y=matrix(1,dim(x)[1],1)){
  return(sum(apply(x,1,mean)*apply(y,1,mean)))
}


compute.SCC.Lemoine <- function(model,
                                damage.type = damage.type,
                                X0,H=500,nb.replic=100,
                                s.values = NaN,
                                coef.multip.values = 1,
                                seed0 = 123,
                                shock = 1){
  # s.values and alpha.values have to be lists
  
  if(is.na(s.values[1])){
    s.values <- list(model$s)
  }
  
  # if(is.na(coef.multip.values[1])){
  #   coef.multip.values <- list(model$damage$coef.multip)
  # }
  
  eta  <- model$eta
  beta <- model$beta
  C0 <- X0[1]
  L0 <- X0[2]
  
  matrix.beta.h <- matrix(exp(-(1-beta)*(1:H)),H,nb.replic*
                            length(s.values)*length(coef.multip.values))
  
  seed <- seed0
  
  all.sim.L <- NULL
  all.sim.T <- NULL
  all.sim.C <- NULL
  all.sim.C.perturb <- NULL
  
  for(s in s.values){
    for(coef.multip in coef.multip.values){
      
      seed <- seed0
      
      model.new   <- model
      model.new$s <- s
      
      set.seed(seed)
      res.simul <- simul.model.Lemoine(model.new,
                                       damage.type = damage.type,
                                       X0,
                                       H=H,nb.replic = nb.replic,
                                       coef.multip = coef.multip)
      set.seed(seed)
      X1 <- X0
      X1[4] <- X0[4] - shock * model$psi1
      X1[5] <- X0[5] - shock * (1 - model$psi1) * model$psi2
      res.simul.perturb <- simul.model.Lemoine(model.new,
                                               damage.type = damage.type,
                                               X0=X1,
                                               H=H,nb.replic = nb.replic,
                                               coef.multip = coef.multip)
      
      all.sim.C         <- cbind(all.sim.C,res.simul$all.C)
      all.sim.C.perturb <- cbind(all.sim.C.perturb,res.simul.perturb$all.C)
      all.sim.L         <- cbind(all.sim.L,res.simul$all.L)
      all.sim.T         <- cbind(all.sim.T,res.simul$all.T)
    }
  }
  
  C.L <- all.sim.C/all.sim.L
  U   <- matrix.beta.h * C.L^(1-eta)/(1-eta)
  U   <- sumE(U)
  
  C.L.perturb <- all.sim.C.perturb/all.sim.L
  U.perturb   <- matrix.beta.h * C.L.perturb^(1-eta)/(1-eta)
  U.perturb   <- sumE(U.perturb)
  
  SCC <- (U.perturb - U)/shock/((C0/L0)^(-eta)/L0) * 10^12/10^9 / 3.667
  
  B         <- (C.L.perturb - C.L)/shock
  M         <- matrix.beta.h * (C.L)^(-eta)/(C0/L0)^(-eta)
  M.perturb <- matrix.beta.h * (C.L.perturb)^(-eta)/(C0/L0)^(-eta)
  MM <- .5*(M + M.perturb)
  SCC.check <- L0 * sumE(MM*B) * 10^12/10^9 / 3.667
  
  NPV <- L0 * sumE(MM,B) * 10^12/10^9 / 3.667
  NPV <- L0 * sum(apply(MM,1,mean) * apply(B,1,mean)) * 10^12/10^9 / 3.667
  
  # Compute Temperature risk premium:
  E.T <- apply(all.sim.T,1,mean)
  E.T.risk.adj <- apply(MM*all.sim.T,1,mean)/apply(MM,1,mean)
  
  return(list(
    U = U,
    U.perturb = U.perturb,
    SCC = SCC,
    SCC.check = SCC.check,
    NPV = NPV,
    MM = MM,
    E.T = E.T,
    E.T.risk.adj = E.T.risk.adj,
    all.sim.C = all.sim.C,
    all.sim.C.perturb = all.sim.C.perturb,
    all.sim.L = all.sim.L,
    all.sim.T = all.sim.T,
    aux.M = apply(MM,1,mean),
    aux.B = apply(B,1,mean)
  ))
}

compute_alternative_damage <- function(T, type,
                                       time_step = 1, # number of years per period
                                       coef.multip = 1){
  if(type=="G. Lemoine"){
    alpha = c(.0018)
    T0  <- .8794
    D <- coef.multip * alpha[1] * time_step * (T - T0)
  }
  if(type=="G. Nordhaus-Sztorc"){
    alpha = .00026
    D <- coef.multip * alpha[1] * time_step * T
  }
  if(type=="G. Dell-Jones-Olken"){
    alpha = .00137
    D <- coef.multip * alpha[1] * time_step * T
  }
  if(type=="G. Weitzman"){
    alpha = c(.000075,3.25)
    D <- coef.multip * alpha[1] * time_step * T^alpha[2]
  }
  if(type=="G. Bansal-Kiku-Ochoa"){
    alpha = c(0.01,0.01,0.05)
    ell0 <- .01
    ell1 <- .01
    d    <- .05
    D <- coef.multip * d * (ell0 + ell1 * T)
  }
  if(type=="L. Nordhaus-Sztorc"){
    alpha = .00266
    D <- coef.multip * alpha*T^2
    D <- 1/(1 + D)
  }
  if(type=="L. Barrage-Nordhaus"){
    alpha = 0.003467
    D <- coef.multip * alpha*T^2
    D <- 1 - D
  }
  if(type=="L. Howard-Sterner"){
    alpha = 0.007438
    D <- coef.multip * alpha*T^2
    D <- 1 - D
  }
  if(type=="L. Weitzman"){
    # Using the Hambel et al. (2021) formulation:
    alpha = c(20.46,6.081,6.754)
    exponent <- alpha[3]
    D <- coef.multip * ((T/alpha[1])^2 + 
                          (T/alpha[2])^exponent)
    D <- 1/(1 + D)
  }
  if(type=="L. Barnett-Brock-Hansen"){
    y.bar <- 2
    alpha = c(0.0127,0.0005,.0005)
    D <- coef.multip * (alpha[1]*T + 
                          alpha[2]*T^2 + 
                          alpha[3] * (T - y.bar)^2 * (T > y.bar))
    D <- exp(- D)
  }
  if(type=="L. Traeger"){
    alpha = c(0.022,1/4)
    xi0 <- alpha[1]
    xi1 <- alpha[2]
    D <- coef.multip * (
      1 - exp(xi0 * (1 - exp(xi1*T)))
    )
    D <- 1 - D
  }
  if(type=="L. Dietz-Stern"){
    alpha = c(0.00284,5.07e-6,6.754)
    pi2      <- alpha[1]
    pi3      <- alpha[2]
    exponent <- alpha[3]
    D <- coef.multip * (pi2*T^2 + pi3*T^exponent)
    D <- 1/(1 + D)
  }
  return(D)
}



