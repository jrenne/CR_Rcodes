
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

Psi <- function(u,v){
  
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
  
  # a <- ell0 * (exp(v*d) - 1) + (u + ell1 * chi * (exp(v*d) - 1)) * Theta * mu +
  #   .5 * (u + ell1 * chi * (exp(v*d) - 1))^2 * (Theta^2*sigma_eta^2 + sigma_zeta^2)
  # b <- (u + ell1 * chi * (exp(v*d) - 1)) * nu
  
  a <- v*mu + u*chi*Theta*mu + (u*chi*Theta + v)^2/2*sigma_eta^2 +
    (u*chi)^2/2*sigma_zeta^2 + ell0*(exp(v*d) - 1)
  b <- u*nu + ell1*(exp(v*d) - 1)
  
  return(list(a=a,
              b=b))
}

phi <- function(x){exp(x)-1}



solve_model <- function(model){
  
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
  
  Ebar <- Theta * mu / (1 - nu)
  Tbar <- Ebar * chi
  
  theta <- (1 - gamma) / (1 - 1/psi)
  
  zbar <- 1
  for(i in 1:50){
    kappa1 <- exp(zbar)/(exp(zbar)-1)
    kappa0 <- kappa1 * zbar - log(exp(zbar)-1)
    
    A1 <- ell1/theta * phi((1-gamma)*d)/(kappa1 - nu)
    A0 <- 1/(kappa1 - 1)*(
      log(delta) + kappa0 + (1 - 1/psi + chi*Theta*A1)*mu + ell0/theta*phi((1-gamma)*d) +
        .5 * theta * ((1 - 1/psi + chi*Theta*A1)^2*sigma_eta^2 + (chi*A1*sigma_zeta)^2)
    )
    zbar <- A0 + A1 * Tbar
  }
  
  ab <- Psi((theta - 1)*A1, - gamma)
  a <- ab$a
  b <- ab$b
  
  rfbar <- - theta * log(delta) - (theta - 1)*kappa0 - 
    (theta - 1)*(1 - kappa1)*A0 - a +
    ((theta - 1)*kappa1*A1 - b) * Tbar
  
  a_m <- theta*log(delta) + (theta - 1)*(kappa0 + (1 - kappa1)*A0)
  b_m <- (theta - 1)*(-kappa1*A1)
  c_m <- (theta - 1)*A1
  d_m <- - gamma
  
  return(list(kappa0 = kappa0,
              kappa1 = kappa1,
              A0 = A0, A1 = A1,
              zbar = zbar,
              rfbar = rfbar,
              Ebar = Ebar, Tbar = Tbar,
              a_m = a_m, b_m = b_m, c_m = c_m, d_m = d_m))
}

compute_ab <- function(model,H){
  a <- 0
  b <- 0
  all_a <- NULL
  all_b <- NULL
  
  model_solution <- solve_model(model)
  a_m <- model_solution$a_m
  b_m <- model_solution$b_m
  c_m <- model_solution$c_m
  d_m <- model_solution$d_m
  
  for(h in 1:H){
    ab <- Psi(c_m+b,d_m)
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


model_solution <- solve_model(model)
Tbar <- model_solution$Tbar

res_prices <- compute_ab(model,10)

plot(res_prices$all_r_a + res_prices$all_r_b*Tbar,type="l")


