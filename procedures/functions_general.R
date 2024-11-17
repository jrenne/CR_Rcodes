
#* CALIBRATION -----------------------------------------------------------------
solveParam4c <- function(model,indic_CRRA=FALSE){
  gamma   <- model$parameters$gamma
  delta   <- model$parameters$delta
  delta_K <- model$parameters$delta_K
  
  if(indic_CRRA==FALSE){# In this case, Epstein-Zin preferences
    model$parameters$A_bar <-
      exp(model$target_vector[["mu_c0"]])/delta - 1 + delta_K
    model$parameters$sigma_a <- model$target_vector[["sigma_c0"]] *
      (model$parameters$A_bar + 1 - delta_K)
  }else{# The CRRA case
    model$parameters$A_bar <-
      exp(model$target_vector[["mu_c0"]]*gamma)/delta - 1 + delta_K
    model$parameters$sigma_a <- model$target_vector[["sigma_c0"]]*gamma *
      (model$parameters$A_bar + 1 - delta_K)
  }
  return(model)
}

solveParam4D <-function(model,alpha=.04){
  T.2100 <- c(2,4)
  tstar <- model$horiz.2100
  Tat0   <- model$vector.ini$ini_Tat
  
  cst     <-log(model$target_vector["stdCumD4"]^2+model$target_vector["ECumD4"]^2)/
    log(model$target_vector["ECumD4"])
  
  mu_D    <-(2-cst)/(2*cst-2)
  
  a_alpha_star <- (tstar*(1 - exp(-alpha)) - (1- exp(-alpha*tstar)))/
    ((1- exp(-alpha*tstar))*(1- exp(-alpha)))
  a_alpha_0 <- ((1 - exp(-alpha*tstar)) - tstar*exp(-alpha*tstar)*(1-exp(-alpha)))/
    ((1- exp(-alpha*tstar))*(1- exp(-alpha)))
  
  ## Linear conditional trajectory:
  # ell.D<-solve(matrix(c(tstar,tstar,
  #                       .5*((tstar-1)*T.2100[1]+(tstar+1)*model$vector.ini$ini_Tat),
  #                       .5*((tstar-1)*T.2100[2]+(tstar+1)*model$vector.ini$ini_Tat)),
  #                     2,2))%*%
  #   matrix((-1-mu_D)/(mu_D)*
  #            c(log(model$target_vector["ECumD2"]),
  #              log(model$target_vector["ECumD4"])),
  #          ncol=1)
  
  ab.D<-solve(matrix(c(tstar,tstar,
                       a_alpha_star*T.2100[1] + a_alpha_0*Tat0,
                       a_alpha_star*T.2100[2] + a_alpha_0*Tat0),2,2))%*%
    matrix((-1-mu_D)*
             c(log(model$target_vector["ECumD2"]),
               log(model$target_vector["ECumD4"])),
           ncol=1)
  
  model$parameters$mu_D <- mu_D
  model$parameters$a_D  <- ab.D[1]
  model$parameters$b_D  <- ab.D[2]
  
  return(model)
}

solveParam4H <-function(model,
                        alpha=.04){
  T.2100 <- c(2,4)
  tstar <- model$horiz.2100
  Tat0   <- model$vector.ini$ini_Tat
  
  C<-matrix(c(model$target_vector["EH2"]-model$vector.ini$ini_H,
              model$target_vector["EH4"]-model$vector.ini$ini_H),
            ncol=1)
  
  a_alpha_star <- (tstar*(1 - exp(-alpha)) - (1- exp(-alpha*tstar)))/
    ((1- exp(-alpha*tstar))*(1- exp(-alpha)))
  a_alpha_0 <- ((1 - exp(-alpha*tstar)) - tstar*exp(-alpha*tstar)*(1-exp(-alpha)))/
    ((1- exp(-alpha*tstar))*(1- exp(-alpha)))
  S<-matrix(c(tstar,tstar,
              a_alpha_star*T.2100[1] + a_alpha_0*Tat0,
              a_alpha_star*T.2100[2] + a_alpha_0*Tat0),2,2)
  
  ## Linear conditional trajectory:
  # S<-matrix(c(.5*((tstar-1)*T.2100[1]+(tstar+1)*model$vector.ini$ini_Tat),
  #             .5*((tstar-1)*T.2100[2]+(tstar+1)*model$vector.ini$ini_Tat),
  #             tstar,tstar),
  #           2,2)
  
  SAT <- solve(S)%*%C
  
  model$parameters$a_H <- SAT[1]
  model$parameters$b_H <- SAT[2]
  
  model$parameters$mu_H  <- model$target_vector["stdH4"]^2/
    (2*model$target_vector["EH4"])
  
  return(model)
}

solveParam4N <-function(model,
                        T.2100 = c(2,4),
                        T.2500 = 4,
                        alpha = .04){
  
  tstar <- model$horiz.2100
  Tat0   <- model$vector.ini$ini_Tat
  
  mu_N   <- model$target_vector["stdCumN4"]^2/(2*model$target_vector["ECumN4"])
  
  kapN  <- matrix(seq(0.6,0.98,length=5000),ncol=1)
  
  C     <- matrix(c(model$target_vector["ECumN2"],
                    model$target_vector["ECumN4"]),
                  ncol=1)
  
  # # Linear parametric trajectory:
  # ellN.condkapN<-apply(kapN,1,function(kN){
  #   N<-matrix(c(rep((1-kN^tstar)/(1-kN),2),
  #               (1-kN^tstar)/(1-kN)*model$vector.ini$ini_Tat+
  #                 (T.2100[1]-model$vector.ini$ini_Tat)/tstar*
  #                 (kN+((tstar-1)*kN-tstar)*kN^tstar)/(1-kN)^2,
  #               (1-kN^tstar)/(1-kN)*model$vector.ini$ini_Tat+
  #                 (T.2100[2]-model$vector.ini$ini_Tat)/tstar*
  #                 (kN+((tstar-1)*kN-tstar)*kN^tstar)/(1-kN)^2),
  #             2,2)*mu_N
  #   ellN  <- matrix(c(0,0),ncol=1)
  #   ellN  <- solve(N)%*%C
  #   return(ellN)})
  
  abN.condkapN<-apply(kapN,1,function(kN){
    b_alpha_star <- 1/(1 - exp(-alpha*tstar))*(
      (1 - kN^tstar)/(1 - kN) -
        (1 - kN^tstar*exp(-alpha*tstar))/(1 - kN*exp(-alpha))
    )
    b_alpha_0 <- - 1/(1 - exp(-alpha*tstar))*(
      (1 - kN^tstar)*exp(-alpha*tstar)/(1 - kN) -
        (1 - kN^tstar*exp(-alpha*tstar))/(1 - kN*exp(-alpha))
    )
    N<-matrix(c(rep((1-kN^tstar)/(1-kN),2),
                b_alpha_0*Tat0 + b_alpha_star*T.2100[1],
                b_alpha_0*Tat0 + b_alpha_star*T.2100[2]),
              2,2)
    
    abN  <- matrix(c(0,0),ncol=1)
    abN  <- solve(N)%*%C
    return(abN)})
  
  # CumN.inf    <- c(1/(1-kapN)) *
  #   c(matrix(c(1,T.2500),nrow=1) %*% abN.condkapN)
  
  # ==========================================
  # ==========================================
  #alpha <- 100
  CumN.inf <- c(apply(rbind(c(1/(1-kapN)),
                            c(T.2500/(1-kapN) - (T.2500 - Tat0)/(1-kapN*exp(-alpha)))
                            #c(T.2500/(1-kapN)) # Previous methodology
  ) * abN.condkapN,2,sum))
  # ==========================================
  # ==========================================
  
  param.CumN  <-rbind(t(kapN),abN.condkapN,t(CumN.inf),
                      rep(model$target_vector["ECumNinf"],length(kapN)),
                      t(CumN.inf-model$target_vector["ECumNinf"]))
  rownames(param.CumN)<-c("kappaN","a.N","b.N","CumN.inf",
                          "Target",
                          "Diff")
  
  abN.kapN  <-param.CumN[,which.min(abs(param.CumN["Diff",]))]
  
  model$parameters$mu_N    <- mu_N
  model$parameters$kappa_N <- abN.kapN["kappaN"]
  model$parameters$a_N     <- abN.kapN["a.N"]
  model$parameters$b_N     <- abN.kapN["b.N"]
  
  return(model)
}


solveParam4T <-function(model){
  # determine mu_T so as to match stdTat2100
  
  print("** Calibration of mu_T by bisection method **")
  
  mu_T.low  <- .00000001
  mu_T.high <- .3
  
  model.low  <- model
  model.high <- model
  model.low$parameters$mu_T  <- mu_T.low
  model.high$parameters$mu_T <- mu_T.high
  model_sol.low  <- model_solve(model.low)
  model_sol.high <- model_solve(model.high)
  
  EV_low  <- EV.fct(model_sol.low,model$horiz.2100)
  EV_high <- EV.fct(model_sol.high,model$horiz.2100)
  
  std.T_at.low  <- sqrt(EV_low$VX$T_at[model$horiz.2100])
  std.T_at.high <- sqrt(EV_high$VX$T_at[model$horiz.2100])
  
  # To be sure the higher bound is high enough:
  while(std.T_at.high<model$target_vector["stdTat2100"]){
    model.high <- model
    model.high$parameters$mu_T <- model.high$parameters$mu_T + .2
    model_sol.high <- model_solve(model.high)
    EV_high <- EV.fct(model_sol.high,model$horiz.2100)
    std.T_at.high <- sqrt(EV_high$VX$T_at[model$horiz.2100])
  }
  
  DIFF     <- abs(model$target_vector["stdTat2100"] - std.T_at.low)
  rel.diff <- DIFF/(std.T_at.high - std.T_at.low)
  
  tol  <- .001
  model.interm  <- model
  while(DIFF > tol){
    
    print(paste(" deviation between higher and lower bound: ",toString(std.T_at.high-std.T_at.low),sep=""))
    
    mu_T.interm <- mu_T.low * .5 + mu_T.high * .5
    
    model.interm$parameters$mu_T  <- mu_T.interm
    model_sol.interm  <- model_solve(model.interm)
    EV_interm         <- EV.fct(model_sol.interm,model$horiz.2100)
    std.T_at.interm   <- sqrt(EV_interm$VX$T_at[model$horiz.2100])
    
    if(std.T_at.interm < model$target_vector["stdTat2100"]){
      model_sol.low <- model.interm
      mu_T.low      <- mu_T.interm
      std.T_at.low  <- std.T_at.interm
    }else{
      model_sol.high <- model.interm
      mu_T.high      <- mu_T.interm
      std.T_at.high  <- std.T_at.interm
    }
    
    DIFF     <- abs(model$target_vector["stdTat2100"] - std.T_at.low)
    rel.diff <- DIFF/(std.T_at.high - std.T_at.low)
  }
  
  return(model_sol.low)
}

#x is a value for the parameter sigma_E
mu_T.bisection <- function(x,model){
  
  model$parameters$mu_T <- x
  model_new <- model_solve(model)
  
  EV_new   <- EV.fct(model_new,model$horiz.2100)
  Var.T_at <- EV_new$VX$T_at[model$horiz.2100]
  
  return(sqrt(Var.T_at)-model$target_vector["stdTat2100"])
}



#----------------Plot with confidence area
make_confidence_area <- function(x,y,pdf_xy,p,tol=10^(-2)){
  # Computes the coordinates of a polygon that delineates a confidence
  #     area associated with probability p.
  # The components of x (and of y) must be equally spaced.
  
  h.x <- x[2] - x[1]
  h.y <- y[2] - y[1]
  
  pdf_xy_h2 <- pdf_xy * h.x * h.y # The sum of all entries of pdf_xy_h2 should be close to 1.
  
  # Use bisection to find limit value of pdf_xy_h2 to get probability p within the area:
  z.up <- max(pdf_xy_h2)
  z.lo <- 0
  error <- 1000
  while(error > tol){
    pdf_xy_h2_aux <- pdf_xy_h2
    pdf_xy_h2_aux[pdf_xy_h2<(z.lo+z.up)/2] <- 0
    aux <- sum(pdf_xy_h2_aux) - p
    if(aux>0){
      z.lo <- (z.lo+z.up)/2
    }else{
      z.up <- (z.lo+z.up)/2
    }
    error <- abs(aux)
    #print(error)
  }
  
  M <- pdf_xy_h2_aux
  M[M>0] <- 1 # Then the matrix M is filled only with 0 and 1.
  x.polygon.1 <- NULL
  x.polygon.2 <- NULL
  y.polygon.1 <- NULL
  y.polygon.2 <- NULL
  for(i in 1:length(x)){
    if(sum(M[i,])>0){
      x.polygon.1 <- c(x.polygon.1,x[i])
      x.polygon.2 <- c(x.polygon.2,x[i])
      index.y <- which(M[i,]>0)
      y.polygon.1 <- c(y.polygon.1,y[index.y[1]])
      y.polygon.2 <- c(y.polygon.2,y[index.y[length(index.y)]])
    }
  }
  
  x.polygon <- c(x.polygon.1,rev(x.polygon.2),x.polygon.1[1])
  y.polygon <- c(y.polygon.1,rev(y.polygon.2),y.polygon.1[1])
  
  return(list(
    x.polygon = x.polygon,
    y.polygon = y.polygon
  ))
}


#----------------Determine confidence areas (for shaded areas)
confidence_intervals_across_horizons <- function(all.Probas,
                                                 values.of.variable,
                                                 nb.values.variable = 200,
                                                 vector.of.CI){
  # all.Probas is of dimension K x H where H is the number of considered
  # horizons. its columns are cdf evaluated at values.of.variable
  # nb.values.variable is the number of points used for the finer grid.
  
  scale.variable.values <- seq(values.of.variable[1],
                               tail(values.of.variable,1),
                               length.out = nb.values.variable+1)
  all.pdf <- NULL
  all.cdf <- NULL
  
  H <- dim(all.Probas)[2]
  
  for(t in 1:H){
    # make sure probas are increasing:
    cdf <- pmax(pmin(all.Probas[,t],1),0)
    for(i in 2:length(cdf)){
      if(cdf[i]<cdf[i-1]){
        cdf[i] <- cdf[i-1]
      }
    }
    
    tmp <- splinefun(x=values.of.variable, y=cdf, method="hyman")
    
    fitted.cdf.values <- tmp(scale.variable.values)
    fitted.pdf.values <- diff(fitted.cdf.values)
    
    all.pdf <- cbind(all.pdf,fitted.pdf.values)
    all.cdf <- cbind(all.cdf,fitted.cdf.values)
  }
  
  # Computation of confidence intervals:
  get.index.CI <- function(q,cdf){
    lower.q <- (1-q)/2
    upper.q <- q + (1-q)/2
    # get index for lower q:
    lower.index <- which(cdf>lower.q)[1]
    upper.index <- which(cdf>upper.q)[1]
    return(c(lower.index,upper.index))
  }
  all.CI <- array(NaN,c(2,dim(all.cdf)[2],length(vector.of.CI)))
  for(i in 1:length(vector.of.CI)){
    q <- vector.of.CI[i]
    for(j in 1:dim(all.cdf)[2]){
      aux <- get.index.CI(q,all.cdf[,j])
      all.CI[1,j,i] <- aux[1] # lower bound
      all.CI[2,j,i] <- aux[2] # upper bound
    }
  }
  
  all.CI <- array(scale.variable.values[all.CI],dim(all.CI))
  
  return(list(scale.variable.values = scale.variable.values,
              all.CI = all.CI,
              all.pdf = all.pdf,
              all.cdf = all.cdf))
}





#----------------Function to extract n'th element of the list
extract<-function(list, n){
  sapply(list, `[`, n)
}

#----------------Laplace transform approximation for damages calibration
LT.CumD <- function(model,u,
                    T0=model$vector.ini$ini_Tat,
                    Tstar,
                    tstar,
                    alpha=.04){
  param <- model$parameters
  
  a_alpha_star <- (tstar*(1 - exp(-alpha)) - (1- exp(-alpha*tstar)))/
    ((1- exp(-alpha*tstar))*(1- exp(-alpha)))
  a_alpha_0 <- ((1 - exp(-alpha*tstar)) - tstar*exp(-alpha*tstar)*(1-exp(-alpha)))/
    ((1- exp(-alpha*tstar))*(1- exp(-alpha)))
  
  res <- u/(1-u*param$mu_D) * (tstar*param$a_D +
                                 param$b_D * (a_alpha_star*Tstar + a_alpha_0*T0))
  return(exp(res))
}

#------- Laplace transform approximation for permafrost-related emissions:
LT.CumN <- function(model,u,
                    T0=model$vector.ini$ini_Tat,
                    Tstar,
                    tstar){
  res <- Param.Gamma0.CumN(model,T0,Tstar,tstar)
  psiCumN <- exp(u*res$mu/(1-u*res$mu) * res$lambda)
  return(psiCumN)
}

#------- Laplace transform approximation for SLR:
LT.SLR <- function(model,u,
                   T0=model$vector.ini$ini_Tat,
                   Tstar,
                   tstar,
                   alpha=.04){
  param <- model$parameters
  
  a_alpha_star <- (tstar*(1 - exp(-alpha)) - (1- exp(-alpha*tstar)))/
    ((1- exp(-alpha*tstar))*(1- exp(-alpha)))
  a_alpha_0 <- ((1 - exp(-alpha*tstar)) - tstar*exp(-alpha*tstar)*(1-exp(-alpha)))/
    ((1- exp(-alpha*tstar))*(1- exp(-alpha)))
  
  res <- u*model$vector.ini$ini_H +
    u/(1-u*param$mu_H) * (tstar*param$a_H +
                            param$b_H * (a_alpha_star*Tstar + a_alpha_0*T0))
  
  return(exp(res))
}

Param.Gamma0.CumN <- function(model,
                              T0 = model$vector.ini$ini_Tat,
                              Tstar,
                              tstar,
                              alpha = .04){
  param  <- model$parameters
  kN     <- param$kappa_N
  
  b_alpha_star <- 1/(1 - exp(-alpha*tstar))*(
    (1 - kN^tstar)/(1 - kN) -
      (1 - kN^tstar*exp(-alpha*tstar))/(1 - kN*exp(-alpha))
  )
  b_alpha_0 <- - 1/(1 - exp(-alpha*tstar))*(
    (1 - kN^tstar)*exp(-alpha*tstar)/(1 - kN) -
      (1 - kN^tstar*exp(-alpha*tstar))/(1 - kN*exp(-alpha))
  )
  
  lambda <- (1-kN^tstar)/(1-kN)*param$a_N/param$mu_N +
    (b_alpha_0 * T0 + b_alpha_star * Tstar)*param$b_N/param$mu_N
  
  return(list(lambda=lambda,
              mu=param$mu_N))
}

Param.Gamma0.H <- function(model,
                           T0 = model$vector.ini$ini_Tat,
                           Tstar,
                           tstar,
                           alpha = .04){
  param  <- model$parameters
  
  a_alpha_star <- (tstar*(1 - exp(-alpha)) - (1- exp(-alpha*tstar)))/
    ((1- exp(-alpha*tstar))*(1- exp(-alpha)))
  a_alpha_0 <- ((1 - exp(-alpha*tstar)) - tstar*exp(-alpha*tstar)*(1-exp(-alpha)))/
    ((1- exp(-alpha*tstar))*(1- exp(-alpha)))
  
  lambda <- model$vector.ini$ini_H/param$mu_H + tstar*param$a_H/param$mu_H +
    (a_alpha_0 * T0 + a_alpha_star * Tstar)*param$b_H/param$mu_H
  
  return(list(lambda=lambda,
              mu=param$mu_H))
}



#------- Fourier transform for permafrost-related emissions:
Fourier.psi <- function(model,gamma,x,
                        T0=model$vector.ini$ini_Tat,
                        Tstar,
                        tstar,
                        psi.LT){
  u       <- c(1i*x)
  psiCumD <- matrix(psi.LT(model,u,T0,Tstar,tstar),ncol=1)
  
  dx <- matrix(x-c(0,x[1:length(x)-1]),length(x),1)
  fx <- outer(x,gamma,
              function(r,c)Im(psiCumD[,1]*exp(-1i*r*c))/r)*dx[,1]
  f  <- 1/2 - 1/pi * apply(fx,2,sum)
  
  f.cdf <- pmax(pmin(f,1),0)
  for(i in 2:length(f)){
    if(f.cdf[i]<f.cdf[i-1]){
      f.cdf[i] <- f.cdf[i-1]}}
  
  return(f.cdf)
}

LT.CumN.infinite <- function(model,u,
                             T0 = model$vector.ini$ini_Tat,
                             Tinf = 4,
                             tstar = NaN,
                             alpha = .04){
  param <- model$parameters
  aux <- u/(1-u*param$mu_N) * 
    (param$a_N / (1-param$kappa_N) +
       param$b_N * Tinf / (1 - param$kappa_N) -
       param$b_N * (Tinf - T0) / (1 - param$kappa_N * exp(-alpha)))
  return(exp(aux))
}


#* SCC -------------------------------------------------------------------------

scc.fct <- function(model_sol,h,
                    all=F,
                    C_0=model_sol$parameters$c0, # in trillions of USD
                    X=model_sol$X,t=0){
  # returns SCC in USD per t CO2 (= per tC / 3.667)
  # if all = TRUE, vector of H \in(0, h).
  
  mat <- which(model_sol$names.var.X=="M_at")
  
  mu.u1.c <- abs(extract(model_sol$mu_u1.t1,mat))
  if(h>(model_sol$Tmax-1)){
    mu.u1.c <- c(mu.u1.c,
                 rep(abs(model_sol$inf_matx$mu_u1[mat]),h-model_sol$Tmax+1)) 
  }
  
  omega_c <- matrix(0,model_sol$n.X,1)
  omega_c[which(model_sol$names.var.X=="Cum_dc")] <- 1
  
  if(h==0){
    scc <- - C_0*10^3 * mu_u.t.fct(model_sol)[[2]][which(model_sol$names.var.X=="M_at")]
  } else{
    scc <- C_0*10^3 * multi.lt.fct(model_sol,omega_c,h,X,t)$uX_t.h * mu.u1.c[h+1]
  }
  
  if(all){# returns all scc until h
    scc.all <- - C_0*mu_u.t.fct(model_sol)[[2]][mat]*10^3
    for(i in 1:(h-1)){
      scc.all[i+1] <- multi.lt.fct(model_sol,omega_c,i,X)$uX_t.h * mu.u1.c[i+1]*C_0*10^3
    }
    return(scc.all/(1-model_sol$parameters$delta)/3.667)
  }
  
  return(scc/(1-model_sol$parameters$delta)/3.667)
}


scc.fct.CRRA <- function(model_sol,
                         C_0=model_sol$parameters$c0, # expressed in trillions of USD
                         X=model_sol$X,
                         t=0,
                         H=200){
  # returns SCC in USD per t CO2 (= per tC / 3.667)
  
  mat     <- which(model_sol$names.var.X=="M_at")
  
  gamma <- model_sol$parameters$gamma
  delta <- model_sol$parameters$delta
  scc <- 0
  scc.decomp <- NULL
  for(h in 1:H){
    Uh <- matrix( (1-gamma) * model_sol$mu_c1,
                  model_sol$n.X,h)
    res.lt <- multi.lt.fct.Uh(model_sol,Uh,X,t)
    psi     <- res.lt$uX_t.h
    psi.Mat <- res.lt$psi.1[mat]
    addit.term <- C_0*10^3/(gamma - 1) * delta^h * psi.Mat * psi
    scc.decomp <- c(scc.decomp,addit.term)
    scc <- scc + addit.term
  }
  
  return(list(
    SCC.C = scc,
    SCC.CO2 = scc/3.667,
    scc.decomp = scc.decomp))
}



#* PRICING ---------------------------------------------------------------------

#--------------Yield curve for ZCB
#* Compute yield curve for a ZCB of maturity h
#* returns a (h*2) matrix, 
#* where the first column corresponds to the date,
#* and the second column corresponds to the rate in percent.
compute_yc<-function(model_sol,h,X=model_sol$X,t=0){
  omega_ZCB <- matrix(0,model_sol$n.X,1)
  yds  <-varphi(model_sol,omega.varphi = omega_ZCB,h,X,t)
  yds.r<-yds$r.t
  yds.d<-yds$date
  
  r<-cbind(yds.d,yds.r)
  
  return(r)
}

#--------------Cst maturity for ZCB
#* Compute constant maturity 'h' for ZCB 'nb' times.
#* returns a (h*2) matrix, 
#* where the first column corresponds to the date,
#* and the second column corresponds to the rate in percent.
compute_cst_h<-function(model_sol,h,nb,X=model_sol$X,t=0){
  EX<-cbind(X,t(extract(EV.fct(model_sol,nb)$EX,(t+1):nb)))
  r<-0
  omega_ZCB <- matrix(0,model_sol$n.X,1)
  for(j in t:(nb-1)){
    yield <-varphi(model_sol,omega.varphi = omega_ZCB,h,EX[,j+1],j)
    r[j+1]<-yield$r.t[h]
  }
  date<-seq(model_sol$tstep*t+model_sol$vec_date[1],
            by=model_sol$tstep,
            length=nb)
  
  dr<-cbind(date,r)
  
  return(dr)
}

#--------------Corollary: TIBs-leverage function
#*i is the variable indexed to the bond
#*H is the maturity
#*T0 is the indexed decided by the issuer at date t
#*chi is the leverage factor
TIB <- function(model_sol,chi,T0,H,i,X=model_sol$X,t=0){
  o.ZCB<-matrix(0,model_sol$n.X,1)
  o.t    <- o.ZCB
  o.t[i] <- 1
  P.tib  <- (1-chi*T0) * varphi(model_sol,o.ZCB,H,X,t)[["P.t"]] +
    chi * varphi.tilde(model_sol,o.t,H,X,t)[["P.tilde"]]
  r.tib <- - log(P.tib)/(model_sol$tstep*(1:H))*100
  
  vec_date<-seq(model_sol$vec_date[2],model_sol$tstep,length=H)
  mylist<-list("P.tib"=P.tib,"r.tib"=r.tib,"date"=vec_date)
  return(mylist)
}

#--------------Proposition: payoff on date t+h = exp(t(omega)%*%X)
#*H is the horizon of the bond [t+1:t+99], max H=98
#*omega.varphi =dim(X)*1 matrix
#*return r in percent
varphi<-function(model_sol,
                 omega.varphi,
                 H,
                 X=model_sol$X,t=0){
  # if((t+H)>=(model_sol$Tmax-1)){
  #   P.pi    <-model_sol$pi
  #   inf     <-rep(list(model_sol$inf_matx$pi.inf),t+H-(model_sol$Tmax-1)+1)
  #   P.pi    <-c(P.pi,inf)
  #   
  #   P.eta_1 <-model_sol$eta1
  #   inf     <-rep(list(model_sol$inf_matx$eta1.inf),t+H-(model_sol$Tmax-1)+1)
  #   P.eta_1 <-c(P.eta_1,inf)
  #   
  #   P.eta_0 <-model_sol$eta0
  #   inf     <-matrix(model_sol$inf_matx$eta0.inf,t+H-(model_sol$Tmax-1)+1,1)
  #   P.eta_0 <-rbind(P.eta_0,inf)
  # }else{
  #   P.pi    <-model_sol$pi
  #   P.eta_1 <-model_sol$eta1
  #   P.eta_0 <-model_sol$eta0
  # }
  
  P.pi    <-model_sol$pi
  P.eta_1 <-model_sol$eta1
  P.eta_0 <-model_sol$eta0
  if((t+H)>=(model_sol$Tmax-1)){# maturity beyond Tmax
    inf     <-rep(list(model_sol$inf_matx$pi.inf),t+H-(model_sol$Tmax-1)+1)
    P.pi    <-c(P.pi,inf)
    inf     <-rep(list(model_sol$inf_matx$eta1.inf),t+H-(model_sol$Tmax-1)+1)
    P.eta_1 <-c(P.eta_1,inf)
    inf     <-matrix(model_sol$inf_matx$eta0.inf,t+H-(model_sol$Tmax-1)+1,1)
    P.eta_0 <-rbind(P.eta_0,inf)
  }
  
  #List of all our U.sh
  U.tk   <- list(P.pi[[t+1]]+omega.varphi)
  P.a.pi <- matrix(NaN,H,1)
  for(h in 1:H){
    if((t+h) <= length(model_sol$pi)){
      P.a.pi[h] <- a1.fct(model_sol,P.pi[[t+h]],t+(h-1))
    }else{
      P.a.pi[h] <- a1.fct.inf(model_sol,P.pi[[t+h]]) 
    }
  }
  
  if(H>1){
    for (h in 2:H){
      uk<-matrix(NaN,nrow=model_sol$n.X,ncol=h)
      for(k in 1:(h-1)){
        if((t+k)<(model_sol$Tmax-1)){
          uk[,k] <- - P.eta_1[[t+k+1]] -
            b1.fct(model_sol,P.pi[[t+k+1]],t+k) + P.pi[[t+k]]
        }else{
          uk[,k] <- - P.eta_1[[t+k+1]] -
            b1.fct.inf(model_sol,P.pi[[t+k+1]]) + P.pi[[t+k]]
        }
      }
      uk[,h]   <- P.pi[[t+h]]+omega.varphi
      U.tk[[h]]<-uk
    }
  }
  
  P.psi<-lapply(1:H,function(h){
    if(h==1){
      multi.lt.fct(model_sol,U.tk[[h]],h,X,t)
    }else{
      multi.lt.fct.Uh(model_sol,U.tk[[h]],X,t)
    }
  })
  
  varphi0<-matrix(NaN,H,1)
  varphi1<-list()
  
  for (h in 1:H){
    varphi0[h]  <--sum(P.eta_0[t+(1:h)])-
      sum(P.a.pi[1:h])+
      P.psi[[h]][[1]]
    if(t>(model_sol$Tmax-2)){
      varphi1[[h]]<--P.eta_1[[t+1]]-
        b1.fct.inf(model_sol,P.pi[[t+1]])+
        P.psi[[h]][[2]]
    }else{
      varphi1[[h]]<--P.eta_1[[t+1]]-
        b1.fct(model_sol,P.pi[[t+1]],t)+
        P.psi[[h]][[2]]
    }
  }
  
  varphi<-matrix(unlist(lapply(1:H,function(h)exp(varphi0[h]+
                                                    t(varphi1[[h]])%*%X))),
                 H,1)
  r.t<--log(varphi)/(model_sol$tstep*(1:H))*100
  
  vec_date<-seq(model_sol$tstep*t+model_sol$vec_date[2],
                by=model_sol$tstep,
                length=H)
  
  mylist<-list("varphi0"=varphi0,"varphi1"=varphi1,"P.t"=varphi,"r.t"=r.t,
               "date"=vec_date)
  return(mylist)
}

#--------------Proposition: payoff on date t+h = exp(t(omega)%*%X)*1_{t(a)%*%X<b}
#*a is a vector dim(X)*1 and b is a scalar associated with the options 
#*a= payment associated with some components of X, and b=strike
varphi.hat<-function(model_sol,omega.v.hat,H,x,a,b,X=model_sol$X,t=0){
  dx<-matrix(x-c(0,x[1:length(x)-1]),length(x),1)
  fx<-matrix(NaN,H,length(x))
  
  for (i in 1:length(x)){
    if(i==length(x)*round(i/length(x),1)){
      print(paste("Progress: ",100*i/length(x),"%",sep=""))
    }
    
    fx[,i]<-Im(varphi(model_sol,1i*a*x[i]+omega.v.hat,H,X,t)[[3]]*
                 exp(-1i*b*x[i]))/x[i]*dx[i]
  }
  print("done")
  varphi.hat<-varphi(model_sol,omega.v.hat,H,X,t)[[3]]/2-
    1/pi*matrix(apply(fx,1,sum),H,1)
  
  r.hat<--log(varphi.hat)/(model_sol$tstep*(1:H))*100
  
  vec_date<-seq(model_sol$tstep*t+model_sol$vec_date[2],
                by=model_sol$tstep,
                length=H)
  
  mylist<-list("P.hat"=varphi.hat,"r.hat"=r.hat,"date"=vec_date)
  
  return(mylist)
}

#--------------Corollary: payoff on date t+h = t(omega)%*%X
varphi.tilde<-function(model_sol,omega.v.tilde,H,X=model_sol$X,t=0){
  eps  <-10^(-6)
  o.ZCB<-matrix(0,model_sol$n.X,1)
  
  varphi.tilde <- (varphi(model_sol,eps*omega.v.tilde,H,X,t)[["P.t"]]-
                     varphi(model_sol,o.ZCB,H,X,t)[["P.t"]])/eps
  
  r.tilde <- - log(varphi.tilde)/(model_sol$tstep*(1:H))*100
  
  vec_date<-seq(model_sol$tstep*t+model_sol$vec_date[2],
                by=model_sol$tstep,
                length=H)
  
  mylist<-list("P.tilde"=varphi.tilde,"r.tilde"=r.tilde,"date"=vec_date)
  
  return(mylist)
}

#--------------Corollary: payoff on date t+h = t(omega)%*%X*1_{t(a)%*%X<b}
#*x is a grid that allow the function to do the approximated Fourier transform
varphi.bar<-function(model_sol,omega.v.bar,H,x,a,b,X=model_sol$X,t=0){
  eps       <-10^-5
  varphi.bar<-(varphi.hat(model_sol,eps*omega.v.bar,H,x,a,b,X,t)[[1]]-
                 varphi.hat(model_sol,0,H,x,a,b,X,t)[[1]])/eps
  
  print("done")
  r.bar<--log(varphi.bar)/(model_sol$tstep*(1:H))*100
  
  vec_date<-seq(model_sol$tstep*t+model_sol$vec_date[2],
                by=model_sol$tstep,
                length=H)
  
  mylist<-list("P.bar"=varphi.bar,"r.bar"=r.bar,"date"=vec_date)
  
  return(mylist)
}


#* FOURIER TRANSFORM -----------------------------------------------------------

#------------------------Fourier defined for cdf (u=0)
#*i is the variable in X we are interested in
#*h is the maturity/horizon
#*u is set to 0 to compute the cdf
#*x is a grid for the integral
#*gamma is the set of possible values for the variable i
psi<-function(model_sol,u,h,i,X=model_sol$X,indic.cum=0){
  U    <-matrix(0,model_sol$n.X,length(u))
  U[i,]<-u
  psi  <-multi.lt.fct.N(model_sol,U,h,X,indic.cum=indic.cum)
  return(psi)
}
fourier<-function(model_sol,x,gamma,h,i,X=model_sol$X,u=0,indic.cum=0){
  # works only for i^th variable
  dx<-matrix(x-c(0,x[1:length(x)-1]),length(x),1)
  s1<-psi(model_sol,u+1i*x,h,i,X,indic.cum=indic.cum)
  fx<-outer(x,gamma,function(r,c)Im(s1[,1]*exp(-1i*r*c))/r)*
    dx[,1]
  f <-1/2-1/pi*apply(fx,2,sum)
  return(f)
}


fourier_complete<-function(model_sol,x,
                           omega,a,b,
                           h,
                           X=model_sol$X,indic.cum=0){
  # works for any linear combination (omega'X) of X variables
  dx<-matrix(x-c(0,x[1:length(x)-1]),length(x),1)
  U <- matrix(omega,model_sol$n.X,length(x))
  U <- U + 1i*matrix(a,ncol=1) %*% matrix(x,nrow=1)
  s1<-multi.lt.fct.N(model_sol,U,h,X,indic.cum=indic.cum)
  fx<-outer(x,b,function(r,c)Im(s1[,1]*exp(-1i*r*c))/r)*
    dx[,1]
  psiomega <- multi.lt.fct.N(model_sol,matrix(omega,ncol=1),
                             h,X,indic.cum=indic.cum)
  f <- c(psiomega)/2 - 1/pi*apply(fx,2,sum)
  return(f)
}


#* MULTI LT, SIMPLE + COMPLETE -------------------------------------------------

#----------------------------------Proposition: Multi-horizon LapT of X
#*U is a matrix of dimension "dim(X)*H", with H>1, the maturity
#*return psi0, psi1 and psi_t.H
multi.lt.fct.Uh <- function(model_sol,Uh,
                            X=model_sol$X,t=0){
  H <- dim(Uh)[2]
  U <- list(Uh[,H,drop=F])
  if((t+H)>(model_sol$Tmax-1)){
    a.h <- matrix(c(a1.fct.inf(model_sol,U[[1]]),rep(NaN,H-1)),H,1) 
    i   <- 1
    #as long as t+H-i is larger than 98=max(a1/b1), we stay time-indep
    while((t+H-i)>(model_sol$Tmax-2)&i<H){
      U[[i+1]] <- Uh[,H-i,drop=F] + b1.fct.inf(model_sol,U[[i]])
      a.h[i+1] <- a1.fct.inf(model_sol,U[[i+1]])
      i<-i+1
    }
    if(i<=(H-1)){
      for (k in i:(H-1)){
        U[[k+1]] <- Uh[,H-k,drop=F] + b1.fct(model_sol,U[[k]],t+H-k)
        a.h[k+1] <- a1.fct(model_sol,U[[k+1]],t+H-(k+1))
      } 
    }
  }else{
    a.h <- matrix(c(a1.fct(model_sol,U[[1]],t+H-1),rep(NaN,H-1)),H,1) 
    if(H>1){
      for (k in 1:(H-1)){
        U[[k+1]] <- Uh[,H-k,drop=F] + b1.fct(model_sol,U[[k]],t+H-k)
        a.h[k+1] <- a1.fct(model_sol,U[[k+1]],t+H-(k+1))
      }
    }
  }
  if(t>(model_sol$Tmax-2)){
    uX_t.h <- exp(apply(a.h,2,sum) + t(b1.fct.inf(model_sol,U[[H]]))%*%X)
    psi.0  <- apply(a.h,2,sum)
    psi.1  <- b1.fct.inf(model_sol,U[[H]])
  }else{
    uX_t.h <- exp(apply(a.h,2,sum) + t(b1.fct(model_sol,U[[H]],t))%*%X)
    psi.0  <- apply(a.h,2,sum)
    psi.1  <- b1.fct(model_sol,U[[H]],t) 
  }
  mylist <- list("psi.0"  = psi.0,
                 "psi.1"  = psi.1,
                 "uX_t.h" = uX_t.h)
  return(mylist)
}

# Corollary: Simple Multihorizon LapT ------------------------------------------
#*h corresponds to the maturity, max 99, then time-independent functions
#*U is a column vector of dimension "dim(X)*1"
#*return psi0, psi1 and psi_t.h
multi.lt.fct <- function(model_sol,U,
                         h,X=model_sol$X,t=0){
  
  param <- model_sol$parameters
  
  U.h   <- list(U)
  a.sum <- 0
  #recursive method with h going to time-indep maturity
  if((t+h)>(model_sol$Tmax-1)){
    i<-1
    while((t+h+1-i)>(model_sol$Tmax-1)){
      a.sum      <- a.sum + a1.fct.inf(model_sol,U.h[[i]])
      U.h[[i+1]] <- b1.fct.inf(model_sol,U.h[[i]])
      
      i<-i+1
    }
    if(i<=h){
      for (k in i:h){
        a.sum      <- a.sum + a1.fct(model_sol,U.h[[k]],t+h-k)
        U.h[[k+1]] <- b1.fct(model_sol,U.h[[k]],t+h-k)
      }  
    }
    uX_t.h <- exp(a.sum+t(U.h[[h+1]])%*%X)
  }else{#Or max of time-dep maturity
    for (k in 1:h){
      a.sum      <- a.sum + a1.fct(model_sol,U.h[[k]],t+h-k)
      U.h[[k+1]] <- b1.fct(model_sol,U.h[[k]],t+h-k)
    }
    uX_t.h <- exp(a.sum+t(U.h[[h+1]])%*%X)
  }
  psi0 <- a.sum
  psi1 <- U.h[[h+1]]
  
  mylist <-   list("psi0"   = psi0,
                   "psi1"   = psi1,
                   "uX_t.h" = uX_t.h)
  return(mylist)
}

#----------------------------------Corollary .N: Simple Multihorizon LapT
#*U is a matrix of dimension "dim(X)*N"
#*return psi_t.h
multi.lt.fct.N <- function(model_sol,U,h,
                           X=model_sol$X,t=0,indic.cum=0){
  
  param <- model_sol$parameters
  
  U.h   <- list(U)
  a.sum <- matrix(0,nrow=dim(U)[2])
  if((t+h)>(model_sol$Tmax-1)){
    i <- 1
    while((t+h+1-i)>(model_sol$Tmax-1)){
      a.sum      <- a.sum+a1.fct.inf(model_sol,U.h[[i]]+indic.cum*(i>1)*U)
      U.h[[i+1]] <- b1.fct.inf(model_sol,U.h[[i]]+indic.cum*(i>1)*U)
      
      i <- i + 1
    }
    if(i<=h){
      for (k in i:h){
        a.sum      <- a.sum + a1.fct(model_sol,U.h[[k]]+indic.cum*(i>1)*U,t+h-k)
        U.h[[k+1]] <- b1.fct(model_sol,U.h[[k]]+indic.cum*(i>1)*U,t+h-k)
      } 
    }
    uX_t.h <- exp(a.sum+t(U.h[[h+1]])%*%X)
    
  }else{
    for (k in 1:h){
      a.sum      <- a.sum+a1.fct(model_sol,U.h[[k]]+indic.cum*(k>1)*U,t+h-k)
      U.h[[k+1]] <- b1.fct(model_sol,U.h[[k]]+indic.cum*(k>1)*U,t+h-k)
    }
    uX_t.h <- exp(a.sum+t(U.h[[h+1]])%*%X)
  }
  return(uX_t.h)
}


#* SOLVING SDF -----------------------------------------------------------------

#  Function that returns mu_u1 and mu_u0 @0 ------------------------------------
#*Tend: General case with Tmax=100, ends up with t=2020
#*If want a specific date --> mu_u.t.fct.all
#*return mu_u0.t, mu_u1.t
mu_u.t.fct <- function(model_sol,
                       Tend = model_sol$Tmax){
  # comes back to initial date if Tend = Tmax
  
  param   <- model_sol$parameters
  
  # Load stationary values:
  mu_c0   <- model_sol$mu_c0
  mu_c1   <- model_sol$mu_c1
  mu_u0.t <- model_sol$mu_u0
  mu_u1.t <- model_sol$mu_u1                                                          
  
  for (i in 2:Tend){
    b.t <- b1.fct(model_sol,
                  (1-param$gamma)*(mu_u1.t + mu_c1), model_sol$Tmax-i)
    a.t <- a1.fct(model_sol,
                  (1-param$gamma)*(mu_u1.t + mu_c1), model_sol$Tmax-i)
    mu_u1.t <- param$delta/(1-param$gamma)*b.t
    mu_u0.t <- param$delta*(mu_u0.t+mu_c0) + param$delta/(1-param$gamma)*a.t
  }
  mylist<-list("mu_u0.t" = mu_u0.t,
               "mu_u1.t" = mu_u1.t)
  return(mylist)
}

# List with all mu_u1.t and mu_u0.t: -------------------------------------------
#*from t=1 to t=Tmax-2, inf(i.e., from t=Tmax-1 onwards) not in list
mu_u.t.fct.all <- function(model_sol){
  
  param     <- model_sol$parameters
  mu_u1.inf <- model_sol$mu_u1
  mu_u0.inf <- model_sol$mu_u0
  Tmax      <- model_sol$Tmax
  
  mu_c1 <- model_sol$mu_c1
  mu_c0 <- model_sol$mu_c0
  
  mu_u1.t         <- list()
  mu_u1.t[[Tmax]] <- mu_u1.inf
  mu_u0.t         <- matrix(NaN,Tmax,1)
  mu_u0.t[Tmax]   <- mu_u0.inf
  for (i in 1:(Tmax-1)){
    mu_u1.t[[Tmax-i]] <- param$delta/(1-param$gamma)*
      b1.fct(model_sol,
             (1-param$gamma)*(mu_u1.t[[Tmax-i+1]] + mu_c1),
             Tmax-i-1)
    mu_u0.t[Tmax-i] <- param$delta*(mu_u0.t[Tmax-i+1] + mu_c0) +
      param$delta/(1-param$gamma)*
      a1.fct(model_sol,
             (1-param$gamma)*(mu_u1.t[[Tmax-i+1]] + mu_c1),
             Tmax-i-1)
  }
  #no inf
  mu_u1.t <- mu_u1.t[-Tmax]
  mu_u0.t <- mu_u0.t[-Tmax]
  #no date 0 (e.g., 2020)
  mu_u1  <- mu_u1.t[-1]
  mu_u0  <- mu_u0.t[-1]
  vec_date <- seq(model_sol$vec_date[2],by=model_sol$tstep,length=Tmax-2)       #infinity not counted
  
  mylist<-list("mu_u0.t1" = mu_u0,
               "mu_u1.t1" = mu_u1,
               "date"     = vec_date)
  return(mylist)
}


#* SOLVE MODEL, MEAN AND SIMUL -------------------------------------------------

# Solve model for infinite mitig mu = 1 ----------------------------------------
#*theta is a set of ini cond. for the optim. of mitig rate mu
model_solve <- function(model,
                        theta       = model$theta0,
                        indic_mitig = TRUE,
                        mu.chosen   = rep(model$parameters$mu0,
                                          model$Tmax),
                        mu_altern   = model$Cum_dc,
                        indic_CRRA  = FALSE){
  # Data preparation
  model_sol <- model
  param     <- model$parameters
  tstep     <- model$tstep
  Tmax      <- model$Tmax
  
  # Create indicators of position of variables (Variables in Z): ---------------
  for(i in 1:model$n.Z){
    eval(parse(text = gsub(" "," ",
                           paste("indic.",model$names.var.X[i],
                                 "<- which(model$names.var.X=='",model$names.var.X[i],
                                 "')",sep=""))))}
  
  # Create indicators of position of variables (Variables in W): ---------------
  for(i in (model$n.Z+1):(model$n.Z+model$n.W)){
    eval(parse(text = gsub(" "," ",
                           paste("indic.",model$names.var.X[i],
                                 "<- which(model$names.var.X=='",model$names.var.X[i],
                                 "')-model$n.Z",sep=""))))}
  
  # Prepare "ell" components: --------------------------------------------------
  ell1.D   <- matrix(0,nrow=model_sol$n.Z+model_sol$n.W)
  ell1.D[indic.T_at] <-
    model_sol$parameters$b_D/model_sol$parameters$mu_D
  model_sol[["ell1.D"]] <- ell1.D
  
  ell1.N   <-matrix(0,nrow=model_sol$n.Z+model_sol$n.W)
  ell1.N[indic.T_at] <-
    model_sol$parameters$b_N/model_sol$parameters$mu_N
  model_sol[["ell1.N"]] <- ell1.N
  
  ell1.T <-matrix(0,nrow=model_sol$n.Z+model_sol$n.W)
  ell1.T[indic.Forc]   <- param$xi_1*model_sol$tstep
  ell1.T[indic.T_at] <-
    1 - param$xi_1*model_sol$tstep*(param$tau/param$nu + param$xi_2)
  ell1.T[indic.T_lo]    <- param$xi_1*model_sol$tstep*param$xi_2
  ell1.T                <- ell1.T/param$mu_T
  model_sol[["ell1.T"]] <- ell1.T
  
  ell1.H                <- matrix(0,nrow=model_sol$n.Z+model_sol$n.W)
  ell1.H[indic.T_at] <- param$b_H/param$mu_H
  model_sol[["ell1.H"]] <- ell1.H
  
  # Construct matrices determining dynamics after Tmax ("infinite"): -----------
  omega.star.inf              <- matrix(0,model_sol$n.Z,model_sol$n.W)
  omega.star.inf[c(indic.delc,indic.y_tilde),indic.eta_A] <-
    param$sigma_a/(param$A_bar+1-param$delta_K)
  if(indic_CRRA){
    omega.star.inf[indic.delc,indic.eta_A] <-
      omega.star.inf[indic.delc,indic.eta_A]/param$gamma
  }
  #omega.star.inf[indic.E_ind,indic.eta_E] <- 0*param$sigma_E
  #omega.star.inf[indic.Forc,indic.eta_F]  <- param$sigma_F
  omega.star.inf[indic.Cum_dc,] <- 
    t(mu_altern$muprice_1[(model_sol$n.Z+1):
                            (model_sol$n.Z+model_sol$n.W)])
  omega.star.inf[indic.delc,indic.D]     <- -1 # Infff
  omega.star.inf[indic.delc,indic.D]     <- 0 # Infff
  if(indic_CRRA){
    omega.star.inf[indic.delc,indic.D] <- 
      omega.star.inf[indic.delc,indic.D]#/param$gamma
  }
  #omega.star.inf[indic.Cum_D,indic.D]    <- -1 
  omega.star.inf[indic.E,indic.N]        <-  1 #shock N
  omega.star.inf[indic.T_at,indic.WT_at] <-  1
  omega.star.inf[indic.H,indic.dH]       <-  1
  #omega.star.inf[indic.T_lo,indic.WT_at] <- param$xi_3*model_sol$tstep

  model_sol[["omega.star.inf"]] <- omega.star.inf
  
  A0.star.inf <- diag(model_sol$n.Z)
  #A0.star.inf[indic.T_at,indic.WT_at] <- - 1
  A0.star.inf[indic.delc,indic.H] <- param$b_sk # Infff
  #A0.star.inf[indic.delc,indic.H] <- 0 # Infff
  if(indic_CRRA){
    A0.star.inf[indic.delc,indic.H] <- 
      A0.star.inf[indic.delc,indic.H]#/param$gamma
  }
  #A0.star.inf[indic.E,indic.E_ind]   <- -1
  A0.star.inf[indic.Forc,indic.M_at] <- -param$tau/(log(2)*param$m_pi*param$m0)
  A0.star.inf[indic.Cum_E,indic.E]   <- -1
  A0.star.inf[indic.Cum_dc,]         <- -t(mu_altern$muprice_1[1:model_sol$n.Z])
  A0.star.inf[indic.Cum_dc,indic.Cum_dc] <- 1
  
  varphi      <- matrix(0,3,3)
  varphi[1,1] <- param$varphi_11
  varphi[2,1] <- param$varphi_12
  varphi[1,2] <- param$varphi_21
  varphi[2,2] <- param$varphi_22
  varphi[3,2] <- param$varphi_23
  varphi[2,3] <- param$varphi_32
  varphi[3,3] <- param$varphi_33
  
  A1.star.inf <- matrix(0,model_sol$n.Z,model_sol$n.Z)
  A1.star.inf[indic.delc,indic.H] <- param$b_sk # Infff
  #A1.star.inf[indic.delc,indic.H] <- 0 # Infff
  if(indic_CRRA){
    A1.star.inf[indic.delc,indic.H] <- A1.star.inf[indic.delc,indic.H]#/param$gamma
  }
  A1.star.inf[indic.y_tilde,indic.y_tilde]    <- 1
  #A1.star.inf[indic.E_ind,indic.y_tilde]      <- 0
  A1.star.inf[indic.M_at,indic.E] <- model_sol$tstep/3.666
  # ============================================
  # ============================================
  A1.star.inf[indic.M_at:indic.M_lo,indic.M_at:indic.M_lo]  <- varphi%^%(model_sol$tstep)                           
  # A1.star.inf[indic.M_at:indic.M_lo,indic.M_at:indic.M_lo]  <- 
  #   diag(3) + model_sol$tstep * (varphi - diag(3))                          
  # ============================================
  # ============================================
  A1.star.inf[indic.T_lo,indic.T_at]     <- param$xi_3*model_sol$tstep
  A1.star.inf[indic.T_lo,indic.T_lo]     <- 1 - param$xi_3*model_sol$tstep
  #A1.star.inf[indic.Cum_D,indic.Cum_D]   <- 1
  A1.star.inf[indic.Cum_E,indic.Cum_E]   <- 1
  A1.star.inf[indic.Cum_dc,indic.Cum_dc] <- 1
  A1.star.inf[indic.H,indic.H]           <- 1
  #A1.star.inf[indic.T_at,indic.T_at]     <- 1
  
  omega0.star.inf <- matrix(0,model_sol$n.Z,1)
  omega0.star.inf[indic.delc,1] <-
    log(param$delta)+log(param$A_bar+1-param$delta_K)
  if(indic_CRRA){
    omega0.star.inf[indic.delc,1] <- omega0.star.inf[indic.delc,1]/param$gamma
  }
  omega0.star.inf[indic.E,1]     <- 0
  #omega0.star.inf[indic.E_ind,1] <- 0
  omega0.star.inf[indic.Forc,1]  <-
    param$tau/log(2)*(log(param$m0)-1)+param$phi_1
  if(length(mu_altern$muprice_0)==1){
    omega0.star.inf[indic.Cum_dc,1] <- mu_altern$muprice_0
  }else{
    omega0.star.inf[indic.Cum_dc,1] <- mu_altern$muprice_0[model$Tmax]
  }
  omega0.star.inf[indic.T_at,1]  <- - param$a_T # lower bound for WT_at
  
  # Compute reduced form of the matrices: 
  A1.inf     <-solve(A0.star.inf) %*% A1.star.inf
  omega0.inf <-solve(A0.star.inf) %*% omega0.star.inf
  omega.inf  <-solve(A0.star.inf) %*% omega.star.inf
  
  model_sol[["varphi"]]      <- varphi
  model_sol[["A0.star.inf"]] <- A0.star.inf
  model_sol[["A1.inf"]]      <- A1.inf
  model_sol[["omega0.inf"]]  <- omega0.inf
  model_sol[["omega.inf"]]   <- omega.inf
  
  # Construct (time-varying deterministic) matrices before Tmax: ---------------
  
  f_ex   <- matrix(NaN, nrow=Tmax, 1)
  E_land <- matrix(NaN, nrow=Tmax, 1)
  gsigma <- matrix(NaN, nrow=Tmax, 1)
  sigma  <- matrix(NaN, nrow=Tmax, 1)
  mu     <- matrix(NaN, nrow=Tmax, 1)
  bp     <- matrix(NaN, nrow=Tmax, 1)
  bc     <- matrix(NaN, nrow=Tmax, 1)
  
  #Radiative forcings
  H2100                <- model_sol$horiz.2100
  f_ex                 <- matrix(rep(param$phi_0,Tmax),Tmax,1)                        
  f_ex[1:H2100]        <- f_ex[1:H2100] + 
    (1/H2100)*(param$phi_1-param$phi_0)*((1:H2100)-1)
  f_ex[(H2100+1):Tmax] <- f_ex[(H2100+1):Tmax] + (param$phi_1-param$phi_0)
  
  #Emissions from deforestation
  E_land[1:Tmax]<- param$eps_0*(1-param$rho)**((1:(Tmax))-1)           
  
  #Carbon intensity
  gsigma[1] <- param$gsigma1                                                      
  for(i in 2:Tmax) gsigma[i] <- gsigma[i-1]*(1+param$delsigma)   
  sigma[1] <- param$sigma0                                                         
  for(i in 2:Tmax) sigma[i] <- sigma[i-1] * exp(gsigma[i-1] * tstep)      
  
  #Abatement cost exogenous components
  bp[1:Tmax] <- param$pback * (1-param$gback)**((1:Tmax)-1)              
  bc[1:Tmax] <- bp[1:Tmax] * sigma[1:Tmax] / param$theta2/1000
  
  #Keep all elements st [1] is t for t=0
  model_sol[["f_ex"]]   <- f_ex
  model_sol[["E_land"]] <- E_land
  model_sol[["gsigma"]] <- gsigma
  model_sol[["sigma"]]  <- sigma
  model_sol[["bp"]]     <- bp
  model_sol[["bc"]]     <- bc
  
  # Construct initial states: --------------------------------------------------
  Z <- matrix(NaN,model_sol$n.Z,1)
  Z[indic.delc] <- log(param$delta) + 
    log((1 - bc[1]*param$mu0**(param$theta2))*param$A_bar + 1 - param$delta_K)
  if(indic_CRRA){
    Z[indic.delc] <- Z[indic.delc]/param$gamma
  }
  Z[indic.y_tilde] <- model_sol$vector.ini$ini_tildey
  Z[indic.E]       <- model_sol$vector.ini$ini_E
  #Z[indic.E_ind]   <- model_sol$vector.ini$ini_Eind
  Z[indic.Forc]    <- model_sol$vector.ini$ini_F
  Z[indic.M_at]    <- model_sol$vector.ini$ini_Mat
  Z[indic.M_up]    <- model_sol$vector.ini$ini_Mup
  Z[indic.M_lo]    <- model_sol$vector.ini$ini_Mlo
  Z[indic.T_at]    <- model_sol$vector.ini$ini_Tat
  Z[indic.T_lo]    <- model_sol$vector.ini$ini_Tlo
  Z[indic.Cum_E]   <- model_sol$vector.ini$ini_CumE
  Z[indic.Cum_dc]  <- model_sol$vector.ini$ini_Cumdelc
  Z[indic.H]       <- model_sol$vector.ini$ini_H
  W <- matrix(0,model_sol$n.W,1)
  W[indic.WT_at] <- model_sol$vector.ini$ini_Tat
  X <- rbind(Z,W)
  model_sol[["X"]]   <- X
  model_sol[["n.X"]] <- length(X)
  
  # Construct mu_c: ------------------------------------------------------------
  mu_c1             <- matrix(0,model_sol$n.Z+model_sol$n.W,1)
  mu_c1[indic.delc] <- 1
  mu_c0             <- 0
  model_sol$mu_c0   <- mu_c0
  model_sol$mu_c1   <- mu_c1
  
  # Solve for "infinite" utility (after Tmax): ---------------------------------
  
  if(indic_CRRA == 0){# This is needed only if Epstein-Zin preferences
    
    # Step 1: construct initial value for muu1_ini (fixed point problem):
    betawGN <- cbind(matrix(0,model_sol$n.Z+model_sol$n.W,model_sol$n.Z),          
                     rbind(
                       matrix(0,model_sol$n.Z,model_sol$n.eta),t(param$Phi),
                       matrix(0,model_sol$n.W-model_sol$n.eta,model_sol$n.eta)
                     ),
                     param$mu_D*model_sol$ell1.D,
                     param$mu_N*model_sol$ell1.N,
                     param$mu_T*model_sol$ell1.T,
                     param$mu_H*model_sol$ell1.H)
    betaGN <- cbind(rbind(t(model_sol$A1.inf),
                          matrix(0,model_sol$n.W,model_sol$n.Z)),
                    matrix(0,model_sol$n.Z+model_sol$n.W,model_sol$n.W)) +
      betawGN %*% rbind(matrix(0,model_sol$n.Z,model_sol$n.Z+model_sol$n.W),
                        cbind(t(model_sol$omega.inf),diag(model_sol$n.W)))
    model_sol[["M"]] <- betaGN
    muu1_ini         <- param$delta*
      solve(diag(model_sol$n.Z+model_sol$n.W) - param$delta*betaGN) %*%
      betaGN %*% mu_c1
    
    # Step 2: Solve for fixed-point problem (muu1_ini):
    RES.fixed.point <- Auxiliary.function(model_sol,muu1_ini)
    mu_u1           <- RES.fixed.point$mu_u1
    mu_u0           <- param$delta/(1-param$delta)*
      mu_c0 + param$delta/((1-param$delta)*(1-param$gamma))*
      a1.fct.inf(model_sol,(1-param$gamma)*(mu_u1+mu_c1))
    
    model_sol[["mu_u1"]] <- mu_u1
    model_sol[["mu_u0"]] <- mu_u0
    model_sol[["ite"]]   <- RES.fixed.point$ite
    model_sol[["dev"]]   <- RES.fixed.point$dev
  }
  
  # Determine (optimal) mitigation path: ---------------------------------------
  
  if(indic_mitig){
    opt <- res.optim(model_sol,theta,
                     Tend = Tmax,
                     X = X,
                     indic_CRRA)
    theta.opt <- opt[[1]]
    mu        <- mu.function(model_sol,theta.opt)
    model_sol[["theta.opt"]] <- theta.opt
    model_sol[["u0"]]        <- abs(opt[[2]])
  }else{# User does not want to optimize over mu_t:
    mu <- mu.chosen
  }
  
  model_sol[["mu"]] <- mu
  
  # Determine model matrices for t < Tmax: -------------------------------------
  model_matrix <- mu_dep(model_sol,model_sol$mu,
                         mu_altern = mu_altern,
                         indic_CRRA = indic_CRRA)
  
  #Keep all elements st [1] is t for t=0
  model_sol[["AC"]]      <- model_matrix$AC
  model_sol[["lambda"]]  <- model_matrix$lambda
  
  #remove first period st starting from t=0, first element is t=1
  model_sol[["A1"]]      <- model_matrix$A1[-1]
  model_sol[["omega"]]   <- model_matrix$omega[-1]
  model_sol[["omega0"]]  <- model_matrix$omega0[-1]
  
  
  # Solve for utility and SDF for t < Tmax -------------------------------------
  
  if(indic_CRRA==FALSE){# Epstein-Zin case
    
    # Determine prices of risk:
    P.mu_u.t1  <- mu_u.t.fct.all(model_sol)
    P.mu_u1.t1 <- P.mu_u.t1[[2]]
    P.pi <- lapply(1:(Tmax-2),                                                    
                   function(t){
                     (1 - param$gamma) * P.mu_u1.t1[[t]] - param$gamma * mu_c1})
    
    # Determine short-term interest rate specification:
    P.eta_0  <- matrix(unlist(lapply(1:(length(P.pi)),function(t){
      -log(param$delta)+
        a1.fct(model_sol,(1-param$gamma)*(P.mu_u1.t1[[t]]+mu_c1),t-1)-
        a1.fct(model_sol,P.pi[[t]],t-1)})), length(P.pi),1)
    P.eta_1  <-lapply(1:(length(P.pi)),
                      function(t){b1.fct(model_sol,
                                         (1-param$gamma)*(P.mu_u1.t1[[t]] + mu_c1),
                                         t-1)-
                          b1.fct(model_sol,P.pi[[t]],t-1)}) #t=0 for b.t
    
    # Short-term rate and prices of risk for t >= Tmax:
    inf_matx <- list("mu_u1" = model_sol$mu_u1)
    
    inf_matx[["pi.inf"]]   <- (1-param$gamma)*inf_matx$mu_u1 - param$gamma*mu_c1
    inf_matx[["eta0.inf"]] <- - log(param$delta) +
      a1.fct.inf(model_sol,(1-param$gamma)*(inf_matx$mu_u1+mu_c1)) -
      a1.fct.inf(model_sol,inf_matx$pi.inf)
    inf_matx[["eta1.inf"]] <- 
      b1.fct.inf(model_sol,(1-param$gamma)*(inf_matx$mu_u1+mu_c1)) - 
      b1.fct.inf(model_sol,inf_matx$pi.inf)
    
    model_sol[["mu_u1.t1"]] <- P.mu_u1.t1
    
  }else{# CRRA case
    
    # Determine prices of risk:
    P.pi <- lapply(1:(Tmax-2),                                                    
                   function(t){- param$gamma*mu_c1})
    
    # Determine short-term interest rate specification:
    P.eta_0  <- matrix(unlist(lapply(1:(length(P.pi)),function(t){
      - log(param$delta) + param$gamma*mu_c0 - a1.fct(model_sol,P.pi[[t]],t-1)})),
      length(P.pi),1)
    P.eta_1  <-lapply(1:(length(P.pi)),
                      function(t){-b1.fct(model_sol,P.pi[[t]],t-1)}) #t=0 for b.t
    
    # Short-term rate and prices of risk for t >= Tmax:
    inf_matx <- list()
    
    inf_matx[["pi.inf"]]   <- - param$gamma*mu_c1
    inf_matx[["eta0.inf"]] <- - log(param$delta) + param$gamma*mu_c0 -
      a1.fct.inf(model_sol,inf_matx$pi.inf)
    inf_matx[["eta1.inf"]] <- - b1.fct.inf(model_sol,inf_matx$pi.inf)
  }
  
  model_sol[["pi"]]       <- P.pi
  model_sol[["eta1"]]     <- P.eta_1
  model_sol[["eta0"]]     <- P.eta_0
  model_sol[["inf_matx"]] <- inf_matx
  
  return(model_sol)
}


# Functions mitig mu optim -----------------------------------------------------
mu.function <- function(model_sol,
                        theta,
                        t.ini = 0,
                        max.mu.date.0 = 1){
  # This returns a vector of size Tmax
  
  if(length(theta)==2){# We use the parametric approach
    
    a <- -log(max.mu.date.0) + abs(theta[1])
    b <- abs(theta[2])
    if(t.ini==0){# modify mu_t since initial period
      mu  <- pmin(exp(- a + b*(1:model_sol$Tmax)),1)
    }else{# in that case, we replace only mu_t's for t>=t.ini
      # This case is used to test the sensitivity of results to assumption
      #    according to which agents do not re-optimize
      mu <- model_sol$mu
      mu[(t.ini+1):model_sol$Tmax] <- pmin(exp(- a + b*((t.ini+1):model_sol$Tmax)),1)
    }
    
  }else{# We directly optimize on the mu_t's
    
    if(length(theta)!= model_sol$Tmax - t.ini){# problem
      print("theta should be of length Tmax - t.ini")
      return(NaN)
    }else{
      mu <- model_sol$mu
      mu[(t.ini+1):model_sol$Tmax] <- exp(theta)/(1+exp(theta))
    }
  }
  
  mu[1] <- model_sol$parameters$mu0
  
  return(mu)
}


# Compute the vectorial repr. of the model + exogenous dependent fct(mu) -------
# mu is an entry argument.
# Outputs are of length Tmax.
mu_dep <- function(model_sol,
                   mu,
                   mu_altern = model_sol$Cum_dc,#used, e.g., in pricing applications
                   indic_CRRA = FALSE){
  
  param <- model_sol$parameters
  Tmax  <- model_sol$Tmax
  
  AC     <- matrix(NaN, nrow=Tmax, 1)
  lambda <- matrix(NaN, nrow=Tmax, 1)
  
  #Exogenous equations independent of mu ---------------------------------------
  #rad. forcings
  f_ex   <- model_sol$f_ex                        
  #emissions from deforestation
  E_land <- model_sol$E_land           
  #Carbon intensity
  gsigma <- model_sol$gsigma
  sigma  <- model_sol$sigma
  #Abatement cost exogenous components
  bp     <- model_sol$bp            
  bc     <- model_sol$bc
  
  #Exogenous equations dependent on mu -----------------------------------------
  #Abatement Cost
  # AC[1]      <- bc[1]     *param$mu0 **(param$theta2) # bc[1]*mu[1]**(param$theta2)#                                   
  # AC[2:Tmax] <- bc[2:Tmax]*mu[2:Tmax]**(param$theta2)
  AC <- bc*mu**(param$theta2)
  
  #lambda
  for (i in 1:Tmax){
    mu_c    <- log(param$delta) + log((1-AC[1:i])*param$A_bar+1-param$delta_K)
    sigma_c <- (1-AC[1:i])*param$sigma_a/
      ((1-AC[1:i])*param$A_bar+1-param$delta_K)
    if(indic_CRRA){
      mu_c    <- mu_c/param$gamma
      sigma_c <- sigma_c/param$gamma
    }
    lambda[i] <- (1-mu[i]) * sigma[i] * param$q0 * exp(sum(mu_c + .5*sigma_c^2))
  }
  
  # Create indicators of position of variables (Variables in Z):
  for(i in 1:model_sol$n.Z){
    eval(parse(text = gsub(" "," ",
                           paste("indic.",model_sol$names.var.X[i],
                                 "<- which(model_sol$names.var.X=='",model_sol$names.var.X[i],
                                 "')",sep=""))))}
  
  # Create indicators of position of variables (Variables in W):
  for(i in (model_sol$n.Z+1):(model_sol$n.Z+model_sol$n.W)){
    eval(parse(text = gsub(" "," ",
                           paste("indic.",model_sol$names.var.X[i],
                                 "<- which(model_sol$names.var.X=='",model_sol$names.var.X[i],
                                 "')-model_sol$n.Z",sep=""))))}
  
  A0.star <- model_sol$A0.star.inf
  #A0.star[indic.T_at,indic.WT_at] <- - 1
  A0.star[indic.delc,indic.H] <- param$b_sk # Infff
  #A0.star[indic.delc,indic.H] <- 0 # Infff
  if(indic_CRRA){
    A0.star[indic.delc,indic.H] <- A0.star[indic.delc,indic.H]#/param$gamma
  }
  
  A1.star<-list()
  for (i in 1:Tmax){
    A1_i <- matrix(0,model_sol$n.Z,model_sol$n.Z)
    A1_i[indic.y_tilde,indic.y_tilde] <- 1
    #A1_i[indic.E_ind,indic.y_tilde]<- lambda[i]
    A1_i[indic.E,indic.y_tilde] <- lambda[i]
    A1_i[indic.M_at,indic.E] <- model_sol$tstep/3.666
    # ============================================
    # ============================================
    A1_i[indic.M_at:indic.M_lo,indic.M_at:indic.M_lo]  <- model_sol$varphi%^%(model_sol$tstep)
    # A1_i[indic.M_at:indic.M_lo,indic.M_at:indic.M_lo]  <- 
    #   diag(3) + model_sol$tstep * (model_sol$varphi - diag(3))                          
    # ============================================
    # ============================================
    A1_i[indic.T_lo,indic.T_at]    <- param$xi_3*model_sol$tstep
    A1_i[indic.T_lo,indic.T_lo]     <- 1-param$xi_3*model_sol$tstep
    #A1_i[indic.Cum_D,indic.Cum_D]  <- 1
    A1_i[indic.Cum_E,indic.Cum_E]   <- 1
    A1_i[indic.Cum_dc,indic.Cum_dc] <- 1
    A1_i[indic.H,indic.H]           <- 1
    #A1_i[indic.T_at,indic.T_at]     <- 1
    A1_i[indic.delc,indic.H]        <-  param$b_sk
    if(indic_CRRA){
      A1_i[indic.delc,indic.H] <- A1_i[indic.delc,indic.H]#/param$gamma
    }
    A1.star[[i]] <- A1_i 
  }
  
  omega0.star<-list()
  for(i in 1:Tmax){
    omega0_i        <- matrix(0,model_sol$n.Z,1)
    if(sum(model_sol$mu_c)==0){
      omega0_i[indic.delc,1]   <- log(param$delta) +
        log((1-AC[i])*param$A_bar + 1 - param$delta_K)
      if(indic_CRRA){
        omega0_i[indic.delc,1] <- omega0_i[indic.delc,1]/param$gamma
      }
    }else{# in that case, we want to impose a growth trajectory
      omega0_i[indic.delc,1]   <- model_sol$mu_c[i]
    }
    omega0_i[indic.E,1]     <- E_land[i] + lambda[i]
    #omega0_i[indic.E_ind,1] <- lambda[i]
    omega0_i[indic.Forc,1]  <- param$tau/log(2)*(log(param$m0)-1)+f_ex[i]
    if(length(mu_altern$muprice_0)==1){
      omega0_i[indic.Cum_dc,1]  <- mu_altern$muprice_0
    }else{
      omega0_i[indic.Cum_dc,1]  <- mu_altern$muprice_0[i]
    }
    omega0_i[indic.T_at,1]  <- - param$a_T # lower bound for WT_at
    omega0.star[[i]]<- omega0_i
  }
  
  omega.star<-list()
  for(i in 1:Tmax){
    omega_i <- matrix(0,model_sol$n.Z,model_sol$n.W)
    omega_i[c(indic.delc,indic.y_tilde),indic.eta_A] <- param$sigma_a*
      (1-AC[i])/((1-AC[i])*param$A_bar+1-param$delta_K)
    if(indic_CRRA){
      omega_i[indic.delc,indic.eta_A] <- omega_i[indic.delc,indic.eta_A]/param$gamma
    }
    omega_i[indic.Cum_dc,] <- t(mu_altern$muprice_1[(model_sol$n.Z+1):
                                                      (model_sol$n.Z+model_sol$n.W)])
    
    omega_i[indic.delc,indic.D]     <- -1 #shock D
    if(indic_CRRA){
      omega_i[indic.delc,indic.D] <- omega_i[indic.delc,indic.D]#/param$gamma
    }
    #omega_i[indic.Cum_D,indic.D]    <- -1 #shock D
    omega_i[indic.E,indic.N]        <-  1 #shock N on E
    omega_i[indic.T_at,indic.WT_at] <-  1
    omega_i[indic.H,indic.dH]       <-  1
    #omega_i[indic.T_lo,indic.WT_at] <- param$xi_3*model_sol$tstep

    omega.star[[i]]      <- omega_i
  }
  A1      <-lapply(1:Tmax,function(i)solve(A0.star)%*%
                     A1.star[[i]])
  omega0  <-lapply(1:Tmax,function(i)solve(A0.star)%*%
                     omega0.star[[i]])
  omega   <-lapply(1:Tmax,function(i)solve(A0.star)%*%
                     omega.star[[i]]) 
  
  return(list("A1"=A1,
              "omega0"=omega0,
              "omega"=omega,
              "AC"=AC,
              "lambda"=lambda))
}



#-------------------------------Results of optimization problem of mu
#*linked to utility function
#*return the value of optimization and the a/b (theta1/2)
res.optim <-function(model_sol,
                     theta,
                     Tend = model_sol$Tmax,
                     X    = model_sol$X,
                     indic_CRRA = FALSE){
  # print(utility.optim(theta,
  #                     model_sol = model_sol,
  #                     Tend = Tend,
  #                     X = X,
  #                     indic_CRRA = indic_CRRA))
  res.optim <- optim(par = theta,
                     utility.optim,
                     model_sol = model_sol,
                     Tend = Tend,
                     X = X,
                     indic_CRRA = indic_CRRA,
                     gr = NULL,
                     method="Nelder-Mead",
                     control=list(maxit = model_sol$MAXIT * sqrt(length(theta))),
                     hessian=FALSE)
  mylist<-list("res.optim$par"=res.optim$par,
               "res.optim$value"=res.optim$value)
  return(res.optim)
}

# Construction of u0 for optim of mitig rate mu: -------------------------------
# This function computes (minus) the utility function at a given date, for a state X,
#  it is aimed to be optimized (over theta) to determine the optimal mitigation path
utility.optim <- function(model_sol,
                          theta,
                          Tend = model_sol$Tmax, # if Tend == Tmax, back to initial period
                          X = model_sol$X,
                          indic_CRRA = FALSE){
  
  # Compute emissions control rate:
  mu <- matrix(NaN, nrow = model_sol$Tmax, 1)
  mu[1:model_sol$Tmax] <- mu.function(model_sol,
                                      theta,
                                      t.ini = model_sol$Tmax - Tend)
  
  # Construct resulting matrix-form model:
  model_matrix <- mu_dep(model_sol, mu,
                         indic_CRRA = indic_CRRA)
  omega0       <- model_matrix$omega0
  omega        <- model_matrix$omega
  A1           <- model_matrix$A1
  
  # Remove date 0:
  model_sol[["A1"]]     <- A1[-1]
  model_sol[["omega"]]  <- omega[-1]
  model_sol[["omega0"]] <- omega0[-1]
  
  if(indic_CRRA==FALSE){
    # In that case, use the Epstein-Zin preferences
    u0 <- compute.utility.EZ(model_sol,Tend,X)
  }else{
    # In that case, use power-utility preferences
    u0 <- compute.utility.CRRA(model_sol)
  }
  
  return(-u0)
}

compute.utility.EZ <- function(model_sol,Tend,X){
  # ============================================================================
  if(Tend == model_sol$Tmax){
    mu_u.1 <- mu_u.t.fct(model_sol,Tend)
  }else{
    #mu_u.1 <- mu_u.t.fct(model_sol,Tend + 1)
    mu_u.1 <- mu_u.t.fct(model_sol,Tend)
  }
  # ============================================================================
  if(is.na(mu_u.1[[1]])){
    return(-10000)
  }
  if(Tend == model_sol$Tmax){
    u0 <- log(model_sol$param$c0) + mu_u.1[[1]] + t(mu_u.1[[2]]) %*% X
  }else{
    u0 <- log(model_sol$param$c0) + mu_u.1[[1]] + t(mu_u.1[[2]]) %*% X
  }
  return(u0)
}


compute.utility.CRRA <- function(model_sol,
                                 X = model_sol$X,
                                 t = 0,
                                 H = 40){
  gamma <- model_sol$parameters$gamma
  delta <- model_sol$parameters$delta
  u0 <- 1 # one unit of date-t consumption
  all.uh <- NULL
  for(h in 1:H){
    Uh <- matrix( (1-gamma) * model_sol$mu_c1,
                  model_sol$n.X,h)
    res.lt <- multi.lt.fct.Uh(model_sol,Uh,X,t)
    u0 <- u0 + delta^h * res.lt$uX_t.h/(1-gamma)
    
    all.uh <- c(all.uh,res.lt$uX_t.h)
  }
  return(u0)
}



#* Function for theoretical mean and variance ----------------------------------
#*h is the end date of estimations, \in(1,99)
#*list starts in start_date+tstep
EV.fct<-function(model_sol,h=NaN){
  if(is.na(h)){
    t <- length(model_sol$vec_date)
  }else{
    t <- h
  }
  
  param <-model_sol$parameters
  if(t>(model_sol$Tmax-1)){
    omega0 <-model_sol$omega0
    inf    <-rep(list(model_sol$omega0.inf),t-model_sol$Tmax+1)
    omega0 <-c(omega0,inf)
    
    omega <-model_sol$omega
    inf   <-rep(list(model_sol$omega.inf),t-model_sol$Tmax+1)
    omega <-c(omega,inf)
    
    A1    <-model_sol$A1
    inf   <-rep(list(model_sol$A1.inf),t-model_sol$Tmax+1)
    A1    <-c(A1,inf)
  }else{
    omega0 <-model_sol$omega0
    omega  <-model_sol$omega
    A1     <-model_sol$A1
  }
  
  X     <-model_sol$X
  n.Z   <-model_sol$n.Z
  n.W   <-model_sol$n.W
  n.eta <-model_sol$n.eta
  
  #Proposition 4 + 5: Conditional mean and variance of W -----------------------
  #alpha1.w, beta1.w, alpha2.w, beta2.w                                         
  alpha1.w    <-lapply(1:t,function(i){
    rbind(matrix(0,n.Z+n.eta,1),
          param$a_D,
          if(i<model_sol$Tmax){
            param$a_N*param$kappa_N^i
          } else{
            param$a_N*0
          },
          param$a_T, # T_at
          param$a_H)
  })
  
  beta1.w   <-lapply(1:t,function(i){
    cbind(matrix(0,n.Z+n.W,n.Z),
          rbind(matrix(0,n.Z,n.eta),
                param$Phi,
                matrix(0,n.W-n.eta,n.eta)),
          matrix(0,n.Z+n.W,n.W-n.eta)
    )+
      param$mu_D*rbind(matrix(0,n.Z+n.eta,n.Z+n.W),
                       t(model_sol$ell1.D),
                       matrix(0,3,n.Z+n.W))+
      (i<model_sol$Tmax)*param$kappa_N^i*param$mu_N*rbind(matrix(0,n.Z+n.eta,n.Z+n.W),
                                                          matrix(0,1,n.Z+n.W),
                                                          t(model_sol$ell1.N),
                                                          matrix(0,2,n.Z+n.W))+
      param$mu_T*rbind(matrix(0,n.Z+n.eta,n.Z+n.W),
                       matrix(0,2,n.Z+n.W),
                       t(model_sol$ell1.T),
                       matrix(0,1,n.Z+n.W))+
      param$mu_H*rbind(matrix(0,n.Z+n.eta,n.Z+n.W),
                       matrix(0,3,n.Z+n.W),
                       t(model_sol$ell1.H))
  })
  
  alpha2.w  <-lapply(1:t,function(i){
    c(rep(1,n.eta),
      2*param$mu_D*param$a_D,
      if(i<model_sol$Tmax){
        2*param$mu_N*param$a_N*param$kappa_N^i
      } else{
        2*param$mu_N*param$a_N*0
      },
      2*param$mu_T*param$a_T,
      2*param$mu_H*param$a_H)*diag(n.W) 
  })
  
  beta2.w.D                      <-matrix(0,n.W,n.W)
  beta2.w.D[(n.eta+1),(n.eta+1)] <-2*param$mu_D^2
  
  matN                  <- matrix(0,n.W,n.W)
  matN[n.eta+2,n.eta+2] <- 2*param$mu_N^2
  beta2.w.N  <-list()
  beta2.w.N  <-lapply(1:t,function(i){
    if(i<model_sol$Tmax){
      matN*param$kappa_N^i
    } else{
      matrix(0,n.W,n.W)
    }
  })
  
  beta2.w.T                      <- matrix(0,n.W,n.W)
  beta2.w.T[(n.eta+3),(n.eta+3)] <- 2*param$mu_T^2
  
  beta2.w.H                      <- matrix(0,n.W,n.W)
  beta2.w.H[(n.eta+4),(n.eta+4)] <- 2*param$mu_H^2
  
  alpha2.w<-lapply(1:t,function(i){
    matrix(alpha2.w[[i]],n.W*n.W,1)
  })
  
  beta2.w <-lapply(1:t,function(i){
    matrix(beta2.w.D,n.W*n.W,1)%*%t(model_sol$ell1.D)+
      matrix(beta2.w.N[[i]],n.W*n.W,1)%*%t(model_sol$ell1.N)+
      matrix(beta2.w.T,n.W*n.W,1)%*%t(model_sol$ell1.T)+
      matrix(beta2.w.H,n.W*n.W,1)%*%t(model_sol$ell1.H)})
  
  EV<-list()
  
  #Proposition 6+7: Conditional mean and variance of X
  #alpha1.k1
  alpha1.k1<-lapply(1:t,function(i){
    rbind(omega0[[i]],matrix(0,model_sol$n.W,1))+
      cbind(matrix(0,model_sol$n.X,model_sol$n.Z),
            rbind(omega[[i]],diag(model_sol$n.W)))%*%alpha1.w[[i]]
  })
  
  #beta1.k1
  beta1.k1<-lapply(1:t,function(i){
    cbind(rbind(A1[[i]],matrix(0,model_sol$n.W,model_sol$n.Z)),
          matrix(0,model_sol$n.X,model_sol$n.W))+
      cbind(matrix(0,model_sol$n.X,model_sol$n.Z),
            rbind(omega[[i]],diag(model_sol$n.W)))%*%beta1.w[[i]]
  })
  
  # Conditional mean for all h -------------------------------------------------
  
  #alpha1.tk
  alpha1.tk<-list(alpha1.k1[[1]])
  if(t>1){
    for (i in 2:t){
      alpha1.tk[[i]]<-alpha1.k1[[i]]+beta1.k1[[i]]%*%(alpha1.tk[[i-1]])
    }
  }
  
  #beta1.tk
  beta1.tk<-list(beta1.k1[[1]])
  if(t>1){
    for (i in 2:t){
      beta1.tk[[i]]<-beta1.k1[[i]]%*%(beta1.tk[[i-1]])
    }
  }
  
  EV[["alpha1.k1"]]<- alpha1.k1
  EV[["beta1.k1"]] <- beta1.k1
  EV[["beta1.w"]]  <- beta1.w
  
  #EXh#
  
  EXh<-list()
  for (i in 1:t){
    EXh[[i]]<-alpha1.tk[[i]]+beta1.tk[[i]]%*%X
  }
  
  EV[["alpha1.tk"]] <- alpha1.tk
  EV[["beta1.tk"]]  <- beta1.tk
  
  EX<-lapply(1:model_sol$n.X,function(i){extract(EXh,i)})
  
  names(EX) <- model_sol$names.var.X
  
  EV[["EX"]]<-EX
  
  # Conditional Variance -------------------------------------------------------
  
  #alpha2.k1
  alpha2.k1<-lapply(1:t,function(i){
    ((t(cbind(diag(model_sol$n.Z),matrix(0,model_sol$n.Z,model_sol$n.W)))%x%
        rbind(diag(model_sol$n.Z),matrix(0,model_sol$n.W,model_sol$n.Z)))%*%
       (omega[[i]]%x%omega[[i]])+
       (t(cbind(matrix(0,model_sol$n.W,model_sol$n.Z),diag(model_sol$n.W)))%x%
          rbind(diag(model_sol$n.Z),matrix(0,model_sol$n.W,model_sol$n.Z)))%*%
       (diag(model_sol$n.W)%x%omega[[i]])+
       (t(cbind(diag(model_sol$n.Z),matrix(0,model_sol$n.Z,model_sol$n.W)))%x%
          rbind(matrix(0,model_sol$n.Z,model_sol$n.W),diag(model_sol$n.W)))%*%
       (omega[[i]]%x%diag(model_sol$n.W))+
       (t(cbind(matrix(0,model_sol$n.W,model_sol$n.Z),diag(model_sol$n.W)))%x%
          rbind(matrix(0,model_sol$n.Z,model_sol$n.W),diag(model_sol$n.W)))%*%
       diag(model_sol$n.W*model_sol$n.W))%*%alpha2.w[[i]]
  })
  
  #beta2.k1
  beta2.k1<-lapply(1:t,function(i){
    ((t(cbind(diag(model_sol$n.Z),matrix(0,model_sol$n.Z,model_sol$n.W)))%x%
        rbind(diag(model_sol$n.Z),matrix(0,model_sol$n.W,model_sol$n.Z)))%*%
       (omega[[i]]%x%omega[[i]])+
       (t(cbind(matrix(0,model_sol$n.W,model_sol$n.Z),diag(model_sol$n.W)))%x%
          rbind(diag(model_sol$n.Z),matrix(0,model_sol$n.W,model_sol$n.Z)))%*%
       (diag(model_sol$n.W)%x%omega[[i]])+
       (t(cbind(diag(model_sol$n.Z),matrix(0,model_sol$n.Z,model_sol$n.W)))%x%
          rbind(matrix(0,model_sol$n.Z,model_sol$n.W),diag(model_sol$n.W)))%*%
       (omega[[i]]%x%diag(model_sol$n.W))+
       (t(cbind(matrix(0,model_sol$n.W,model_sol$n.Z),diag(model_sol$n.W)))%x%
          rbind(matrix(0,model_sol$n.Z,model_sol$n.W),diag(model_sol$n.W)))%*%
       diag(model_sol$n.W*model_sol$n.W))%*%beta2.w[[i]]
  })
  
  #alpha2.tk
  alpha2.tk<-list(alpha2.k1[[1]])
  if(t>1){
    for (i in 2:t){
      alpha2.tk[[i]]<-alpha2.k1[[i]]+beta2.k1[[i]]%*%alpha1.tk[[i]]+
        (beta1.k1[[i]]%x%beta1.k1[[i]])%*%alpha2.tk[[i-1]]
    }
  }
  
  #beta2.tk
  beta2.tk<-list(beta2.k1[[1]])
  if(t>1){
    for (i in 2:t){
      beta2.tk[[i]]<-beta2.k1[[i]]%*%beta1.tk[[i]]+
        (beta1.k1[[i]]%x%beta1.k1[[i]])%*%beta2.tk[[i-1]]
    }
  }
  
  vecVXh  <-lapply(1:t,
                   function(x)alpha2.tk[[x]]+beta2.tk[[x]]%*%X)
  VXh     <-lapply(1:t,
                   function(x)matrix(vecVXh[[x]],
                                     model_sol$n.X,model_sol$n.X))
  VXh.diag<-lapply(1:t,
                   function(x)diag(VXh[[x]]))
  VX      <-lapply(1:model_sol$n.X, 
                   function(i)extract(VXh.diag,i))
  
  names(VX) <- model_sol$names.var.X
  
  #lower and upper bound
  upper   <-lapply(1:model_sol$n.X,function(x){EX[[x]]+2*(VX[[x]])^(1/2)})
  lower   <-lapply(1:model_sol$n.X,function(x){EX[[x]]-2*(VX[[x]])^(1/2)})
  
  bounds  <-lapply(1:model_sol$n.X,function(x){rbind(lower[[x]],upper[[x]])})
  
  vec_date<-seq(model_sol$vec_date[2],by=model_sol$tstep,length=t)
  
  EV[["CovX"]]  <-VXh
  EV[["VX"]]    <-VX
  EV[["bounds"]]<-bounds
  EV[["date"]]  <-vec_date
  EV[["EXh"]] <- EXh
  return(EV)
}


#* Simulations -----------------------------------------------------------------
#CHANGE
#------------------Function to simulate state variables of our model
#*New approximation of the radiative forcings (linear)
#*nb.simul is number of periods for 1 simulation \in(1,99)
#*nb.traj is the number of simulations done for nb.simul periods
#*for t>99, time-independent matrices
simul.function<-function(model_sol,nb.simul.t,nb.traj,setseed=NaN){
  
  if(!is.na(setseed)){
    set.seed(setseed)
  }
  
  #useful parameters
  nb.simul<-nb.simul.t+1 #take into account t=0
  
  tstep<-model_sol$tstep
  t    <-1:nb.simul
  param<-model_sol$parameters
  
  n.Z   <- model_sol$n.Z
  n.W   <- model_sol$n.W
  n.X   <- model_sol$n.X
  n.eta <- model_sol$n.eta
  
  if(nb.simul>model_sol$Tmax){
    omega0<-model_sol$omega0
    inf   <-rep(list(model_sol$omega0.inf),nb.simul-model_sol$Tmax)
    omega0<-c(omega0,inf)
    
    omega <-model_sol$omega
    inf   <-rep(list(model_sol$omega.inf),nb.simul-model_sol$Tmax)
    omega <-c(omega,inf)
    
    A1    <-model_sol$A1
    inf   <-rep(list(model_sol$A1.inf),nb.simul-model_sol$Tmax)
    A1    <-c(A1,inf)
  }else{
    omega0<-model_sol$omega0
    omega <-model_sol$omega
    A1    <-model_sol$A1
  }
  
  eta <- matrix(0,n.eta,nb.traj)
  
  # Create indicators of position of variables (Variables in Z):
  for(i in 1:n.Z){
    eval(parse(text = gsub(" "," ",
                           paste("indic.",model_sol$names.var.X[i],
                                 "<- which(model_sol$names.var.X=='",model_sol$names.var.X[i],
                                 "')",sep=""))))}
  
  # Create indicators of position of variables (Variables in W):
  for(i in (n.Z+1):(n.Z+n.W)){
    eval(parse(text = gsub(" "," ",
                           paste("indic.",model_sol$names.var.X[i],
                                 "<- which(model_sol$names.var.X=='",model_sol$names.var.X[i],
                                 "')-n.Z",sep=""))))}
  
  Z <- list(matrix(model_sol$X[1:n.Z],n.Z,nb.traj))
  W <- list(matrix(model_sol$X[(n.Z+1):n.X],n.W,nb.traj))
  X <- list(rbind(Z[[1]],W[[1]]))
  
  for (i in 2:nb.simul) {
    W[[i]] <- matrix(0,n.W,nb.traj)
    
    eta <- param$Phi %*% eta + matrix(rnorm(n.eta * nb.traj),n.eta,nb.traj)
    W[[i]][1:n.eta,] <- eta
    
    W[[i]][indic.D,] <-
      rgamma(nb.traj,rpois(nb.traj,
                           pmax(0,(param$a_D/param$mu_D+
                                     t(model_sol$ell1.D)%*%
                                     rbind(Z[[i-1]],W[[i-1]])))),scale=param$mu_D)
    
    W[[i]][indic.N,] <-
      rgamma(nb.traj,rpois(nb.traj,pmax(0,(param$kappa_N^(i-1)*param$a_N/param$mu_N+
                                             param$kappa_N^(i-1)*t(model_sol$ell1.N)%*%
                                             rbind(Z[[i-1]],W[[i-1]])))),scale=param$mu_N)
    
    W[[i]][indic.WT_at,] <-
      rgamma(nb.traj,rpois(nb.traj,pmax(0,(param$a_T/param$mu_T+
                                             t(model_sol$ell1.T)%*%
                                             rbind(Z[[i-1]],W[[i-1]])))),scale=param$mu_T)
    
    W[[i]][indic.dH,] <-
      rgamma(nb.traj,rpois(nb.traj,pmax(0,(param$a_H/param$mu_H+
                                             t(model_sol$ell1.H)%*%
                                             rbind(Z[[i-1]],W[[i-1]])))),scale=param$mu_H)
    
    # Z variables:
    Z_1.0  <- apply(Z[[i-1]],2,
                    function(x){A1[[i-1]]%*%x+omega0[[i-1]]})
    Wz     <- apply(W[[i]],2,function(x){omega[[i-1]]%*%x})
    Z[[i]] <- Z_1.0+Wz
    X[[i]] <- rbind(Z[[i]],W[[i]])
  }
  
  #extract path over time and remove time 0 (i.e., X_0)
  delc   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[indic.delc,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  y_tilde<-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[indic.y_tilde,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  E      <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[indic.E,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  # E_ind  <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[indic.E_ind,]),1:nb.traj)),
  #                 ncol=nb.traj)[-1,]
  Forc   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[indic.Forc,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  M_at   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[indic.M_at,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  M_up   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[indic.M_up,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  M_lo   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[indic.M_lo,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  T_at   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[indic.T_at,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  T_lo   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[indic.T_lo,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  Cum_E  <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[indic.Cum_E,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  Cum_dc <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[indic.Cum_dc,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  H      <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[indic.H,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  
  D    <- matrix(t(extract(lapply(W[1:nb.simul],function(x) x[indic.D,]),1:nb.traj)),
                 ncol=nb.traj)[-1,]
  N    <- matrix(t(extract(lapply(W[1:nb.simul],function(x) x[indic.N,]),1:nb.traj)),
                 ncol=nb.traj)[-1,]
  WT_at <- matrix(t(extract(lapply(W[1:nb.simul],function(x) x[indic.WT_at,]),1:nb.traj)),
                 ncol=nb.traj)[-1,]
  dH   <- matrix(t(extract(lapply(W[1:nb.simul],function(x) x[indic.dH,]),1:nb.traj)),
                 ncol=nb.traj)[-1,]
  
  X      <- X[-1]
  
  vec_date <- seq(model_sol$vec_date[2],
                  by=model_sol$tstep,length=nb.simul-1)
  
  
  mylist<-c(list("delc"=delc,"y_tilde"=y_tilde,"E"=E,"Forc"=Forc,
                 "M_at"=M_at,"M_up"=M_up,"M_lo"=M_lo,"T_lo"=T_lo,"T_at"=T_at,
                 "Cum_E"=Cum_E,"Cum_dc"=Cum_dc,"H"=H),
            list("D"=D,"N"=N,"WT_at"=WT_at,"dH"=dH),
            list("X"=X))
  
  return(mylist)
}


#* LAPLACE TRANSFORMS (LapT) FCTS ----------------------------------------------
#CHANGE if add gamma0
#-------------------------------------Proposition: LapT of W
#*Functions for 'a.w.t' and 'b.w.t'
#*U is a dim(W)*N vector
#*

a1.w.fct <- function(model_sol,U,t){
  param <-model_sol$parameters
  n.eta <-model_sol$n.eta
  
  u.eta<-matrix(U[1:n.eta,],n.eta,dim(U)[2])
  u.D  <- U[n.eta+1,]
  u.N  <- U[n.eta+2,]
  u.T  <- U[n.eta+3,]
  u.H  <- U[n.eta+4,]
  
  a.w <- 0
  a.w[1:dim(U)[2]] <- 0.5*matrix(apply(u.eta*u.eta,2,sum),
                                 dim(U)[2],1)+
    matrix(u.D*param$a_D/(1-u.D*param$mu_D),
           dim(U)[2],1)+
    #===========================================================================
  #===========================================================================
  #===========================================================================
  # matrix((param$kappa_N^(t+1))*u.N*param$a_N/(1-u.N*param$mu_N),
  #        dim(U)[2],1) + 
  #===========================================================================
  #===========================================================================
  #===========================================================================
  matrix((param$kappa_N^t)*u.N*param$a_N/(1-u.N*param$mu_N),
         dim(U)[2],1) + 
    matrix(u.T*param$a_T/(1-u.T*param$mu_T),
           dim(U)[2],1) + 
    matrix(u.H*param$a_H/(1-u.H*param$mu_H),
           dim(U)[2],1)
  
  a.w <- matrix(a.w,dim(U)[2],1)
  
  if(!is.complex(u.D)|!is.complex(u.N)){
    a.w[u.D >= 1/param$mu_D] <- NaN
    a.w[u.N >= 1/param$mu_N] <- NaN
    a.w[u.T >= 1/param$mu_T] <- NaN
  }
  
  return(a.w)
}

b1.w.fct<-function(model_sol,U,t){
  param <-model_sol$parameters
  n.eta <-model_sol$n.eta
  
  u.eta<-matrix(U[1:n.eta,],n.eta,dim(U)[2])
  u.D  <-matrix(U[n.eta+1,],ncol=dim(U)[2])
  u.N  <-matrix(U[n.eta+2,],ncol=dim(U)[2])
  u.T  <-matrix(U[n.eta+3,],ncol=dim(U)[2])
  u.H  <-matrix(U[n.eta+4,],ncol=dim(U)[2])
  
  a <- matrix(0,
              model_sol$n.Z+model_sol$n.W,dim(U)[2])
  a[(model_sol$n.Z+1):
      (model_sol$n.Z+n.eta),]  <- t(param$Phi)%*%u.eta
  
  b <- matrix(model_sol$ell1.D%o%
                (u.D*param$mu_D/(1-u.D*param$mu_D))+
                model_sol$ell1.N%o%
                # ==============================================================
              # ==============================================================
              # ==============================================================
              # (u.N*param$mu_N*(param$kappa_N^(t+1))/
              #    (1-u.N*param$mu_N))+
              # ==============================================================
              # ==============================================================
              # ==============================================================
              (u.N*param$mu_N*(param$kappa_N^t)/
                 (1-u.N*param$mu_N))+
                model_sol$ell1.T%o%
                (u.T*param$mu_T/(1-u.T*param$mu_T))+
                model_sol$ell1.H%o%
                (u.H*param$mu_H/(1-u.H*param$mu_H)),
              model_sol$n.Z+model_sol$n.W,
              dim(U)[2])
  
  b.w <- a + b
  
  if(!is.complex(u.D)|!is.complex(u.N)){
    b.w[,u.D >= 1/param$mu_D] <- NaN
    b.w[,u.N >= 1/param$mu_N] <- NaN
    b.w[,u.T >= 1/param$mu_T] <- NaN
  }
  
  return(b.w)
}


#----Proposition: One-Period ahead Laplace transform of X ----------------------
#*MIN value for t: 0; MAX value for t: Tmax-2
#*U: n.X*N matrix

a1.fct <- function(model_sol,U,t){
  # t has to be < Tmax - 1 (otherwise use a1.fct.inf)
  omega.0 <- model_sol$omega0
  omega   <- model_sol$omega
  
  Uz <- matrix(U[1:model_sol$n.Z,],model_sol$n.Z,dim(U)[2])
  Uw <- matrix(U[(model_sol$n.Z+1):(model_sol$n.Z+model_sol$n.W),],
               model_sol$n.W,dim(U)[2])
  
  a <- t(Uz)%*%omega.0[[t+1]] +
    a1.w.fct(model_sol,Uw + t(omega[[t+1]])%*%Uz, t)
  
  return(a)
}

b1.fct <- function(model_sol,U,t){
  # t has to be < Tmax - 1 (otherwise use b1.fct.inf)
  A.1   <- model_sol$A1
  omega <- model_sol$omega
  
  Uz <- matrix(U[1:model_sol$n.Z,],model_sol$n.Z,dim(U)[2])
  Uw <- matrix(U[(model_sol$n.Z+1):model_sol$n.X,],model_sol$n.W,dim(U)[2])
  
  b <- rbind(t(A.1[[t+1]]) %*% Uz,matrix(0,model_sol$n.W,dim(U)[2])) +
    b1.w.fct(model_sol,Uw + t(omega[[t+1]])%*%Uz, t)
  
  return(b)
}


#* INFINITE AND GN -------------------------------------------------------------

#---- Functions for infinite a,b
#*U: matrix dim(U)=n.X*N
a1.w.fct.inf<-function(model_sol,U){
  param <-model_sol$parameters
  n.eta <-model_sol$n.eta
  
  u.eta<-matrix(U[1:n.eta,],n.eta,dim(U)[2])
  u.D  <- U[n.eta+1,]
  u.N  <- U[n.eta+2,]
  u.T  <- U[n.eta+3,]
  u.H  <- U[n.eta+4,]
  
  a.w<-0
  a.w[1:dim(U)[2]]<- 0.5*matrix(apply(u.eta*u.eta,2,sum),
                                dim(U)[2],1)+
    matrix(u.D*param$a_D/(1-u.D*param$mu_D),dim(U)[2],1)+
    matrix(u.T*param$a_T/(1-u.T*param$mu_T),dim(U)[2],1)+
    matrix(u.H*param$a_H/(1-u.H*param$mu_H),dim(U)[2],1)
  # Note: - nothing for N (no additional permafrost releases after Tmax)
  #       - nothing for T (a_T = 0)
  
  a.w<-matrix(a.w,dim(U)[2],1)
  
  if(!is.complex(u.D)|!is.complex(u.N)){
    a.w[u.D >= 1/param$mu_D]     <- NaN
    a.w[u.N >= 1/param$mu_N]     <- NaN
    a.w[u.T >= 1/param$mu_T]     <- NaN
  }
  
  return(a.w)
}

b1.w.fct.inf <- function(model_sol,U){
  param <- model_sol$parameters
  n.eta <- model_sol$n.eta
  
  u.eta <- matrix(U[1:n.eta,],n.eta,dim(U)[2])
  u.D   <- matrix(U[n.eta+1,],ncol=dim(U)[2])
  u.N   <- matrix(U[n.eta+2,],ncol=dim(U)[2])
  u.T   <- matrix(U[n.eta+3,],ncol=dim(U)[2])
  u.H   <- matrix(U[n.eta+4,],ncol=dim(U)[2])
  
  a <- matrix(0,
              model_sol$n.Z+model_sol$n.W,dim(U)[2])
  a[(model_sol$n.Z+1):
      (model_sol$n.Z+n.eta),]  <- t(param$Phi)%*%u.eta
  
  b <- matrix(model_sol$ell1.D %o% (u.D*param$mu_D/(1-u.D*param$mu_D)) +
                model_sol$ell1.T %o% (u.T*param$mu_T/(1-u.T*param$mu_T)) +
                model_sol$ell1.H %o% (u.H*param$mu_H/(1-u.H*param$mu_H)),
              model_sol$n.Z+model_sol$n.W,dim(U)[2])
  # Note: - nothing for N (no additional permafrost releases after Tmax)
  
  b.w <- a + b
  
  if(u.D >= 1/param$mu_D){print("Warning D")}
  if(u.N >= 1/param$mu_N){print("Warning N")}
  if(u.T >= 1/param$mu_T){print("Warning T")}
  
  return(b.w)
}

a1.fct.inf <- function(model_sol,U){
  omega.0.inf <- model_sol$omega0.inf
  omega.inf   <- model_sol$omega.inf
  
  Uz <- matrix(U[1:model_sol$n.Z,],model_sol$n.Z,dim(U)[2])
  Uw <- matrix(U[(model_sol$n.Z+1):(model_sol$n.Z+model_sol$n.W),],model_sol$n.W,
               dim(U)[2])
  
  a.inf <- t(Uz)%*%omega.0.inf +
    a1.w.fct.inf(model_sol,Uw+t(omega.inf)%*%Uz)
  
  return(a.inf)
}

b1.fct.inf <- function(model_sol,U){
  omega.inf <- model_sol$omega.inf
  A.1.inf   <- model_sol$A1.inf
  
  Uz <- matrix(U[1:model_sol$n.Z],model_sol$n.Z,1)
  Uw <- matrix(U[(model_sol$n.Z+1):(model_sol$n.Z+model_sol$n.W)],model_sol$n.W,1)
  
  b.inf <- rbind(t(A.1.inf) %*% Uz,matrix(0,model_sol$n.W,1)) +
    b1.w.fct.inf(model_sol,Uw + t(omega.inf) %*% Uz)
  
  return(b.inf)
}


# Solution method for mu_u1 (fixed point) --------------------------------------
#*x0: matrix of dimension "n.X*1", initial guess of infinite values
#*return nb of iterations, dev, mu_u1.inf, list dev

Auxiliary.function <- function(model_sol,x0){
  param     <- model_sol$parameters
  A1.inf    <- model_sol$A1.inf
  omega.inf <- model_sol$omega.inf
  #Jacobian
  mu_u1   <- x0
  eps     <- param$eps.GN
  J       <- matrix(0,model_sol$n.Z+model_sol$n.W,model_sol$n.Z+model_sol$n.W)
  dev     <- matrix(1,model_sol$n.Z+model_sol$n.W,1)
  ite     <- 0
  listdev <- list(dev)
  tol     <- max(abs(dev*param$tol.GN))
  while((max(abs(dev))>tol)&(ite<50)){
    b     <-  b1.fct.inf(model_sol,(1-param$gamma)*(mu_u1 + model_sol$mu_c1))
    dev   <- -mu_u1+param$delta/(1-param$gamma)*b
    mu_u1 <-  param$delta/(1-param$gamma)*b
    ite   <-  ite + 1
    listdev[[ite]] <- dev
  }
  mylist <- list("ite"=ite,"dev"=dev,"mu_u1"=mu_u1,"listdev"=listdev)
  return(mylist)
}



### New pricing functions where varphi is evaluated in parallel ----------------
### at different omega's. For only one maturity at a time.


multi.lt.fct.Uh_MULTI <- function(model_sol,Uh,X=model_sol$X,t=0){
  H <- dim(Uh)[3] # number of horizons
  K <- dim(Uh)[2]
  U <- list(matrix(Uh[,,H],ncol=K))
  if((t+H) > (model_sol$Tmax-1)){
    a.h <- matrix(rbind(c(a1.fct.inf(model_sol,U[[1]])),
                        matrix(NaN,H-1,K)),H,K)   
    i <- 1
    #as long as t+H-i is larger than 98=max(a1/b1), we stay time-indep
    while((t+H-i) > (model_sol$Tmax-2)&i<H){
      U[[i+1]] <- Uh[,,H-i]+b1.fct.inf(model_sol,U[[i]])
      a.h[i+1] <- c(a1.fct.inf(model_sol,U[[i+1]]))
      i<-i+1
    }
    if(i <= (H-1)){
      for (k in i:(H-1)){
        U[[k+1]] <- Uh[,,H-k]+b1.fct(model_sol,U[[k]],t+H-k)
        a.h[k+1] <- c(a1.fct(model_sol,U[[k+1]],t+H-(k+1)))
      } 
    }
  }else{
    a.h <- matrix(rbind(c(a1.fct(model_sol,U[[1]],t+H-1)),
                        matrix(NaN,H-1,K)),H,K) 
    if(H>1){
      for (k in 1:(H-1)){
        U[[k+1]]  <- Uh[,,H-k]+b1.fct(model_sol,U[[k]],t+H-k)
        a.h[k+1,] <- c(a1.fct(model_sol,U[[k+1]],t+H-(k+1)))
      }
    }
  }
  if(t>(model_sol$Tmax-2)){
    uX_t.h <- exp(apply(a.h,2,sum)+t(b1.fct.inf(model_sol,U[[H]]))%*%X)
    psi.0  <- apply(a.h,2,sum)
    psi.1  <- b1.fct.inf(model_sol,U[[H]])
  }else{
    uX_t.h <- exp(apply(a.h,2,sum)+t(b1.fct(model_sol,U[[H]],t))%*%X)
    psi.0  <- apply(a.h,2,sum)
    psi.1  <- b1.fct(model_sol,U[[H]],t) 
  }
  mylist <- list("psi.0"  = psi.0,
                 "psi.1"  = psi.1,
                 "uX_t.h" = uX_t.h)
  return(mylist)
}


varphi_multi_omega<-function(model_sol,omega,H,X=model_sol$X,t=0){
  # omega of dimension n.X x K
  
  K <- dim(omega)[2]
  
  P.pi <- t(matrix(unlist(model_sol$pi),
                   model_sol$n.X,model_sol$Tmax-2))
  eta0 <- model_sol$eta0
  eta1 <- t(matrix(unlist(model_sol$eta1),
                   model_sol$n.X,model_sol$Tmax-2))
  
  P.a.pi <- matrix(NaN,H-1,1)
  P.b.pi <- matrix(NaN,H-1,model_sol$n.X)
  
  U <- array(0,c(model_sol$n.X,K,H))
  
  if(H>1){
    for(k in 1:(H-1)){
      P.a.pi[k]  <- a1.fct(model_sol,matrix(P.pi[t+k+1,],ncol=1),t+k)
      P.b.pi[k,] <- b1.fct(model_sol,matrix(P.pi[t+k+1,],ncol=1),t+k)
      
      U[,,k] <- - eta1[t+k+1,] - P.b.pi[k,] + P.pi[t+k,]
    }
  }
  U[,,H] <- P.pi[t+H,] + omega
  
  res.psi <- multi.lt.fct.Uh_MULTI(model_sol, Uh=U, X=X, t=t)
  
  # Specific treatment for k=1:
  P.a.pi.1 <- a1.fct(model_sol,matrix(P.pi[t+1,],ncol=1),t)
  P.b.pi.1 <- b1.fct(model_sol,matrix(P.pi[t+1,],ncol=1),t)
  
  varphi0 <- - sum(eta0[(t+1):(t+H)]) - c(P.a.pi.1) - sum(P.a.pi) +
    res.psi$psi.0
  varphi1 <- matrix(- eta1[t+1,] - P.b.pi.1,model_sol$n.X,K) + 
    res.psi$psi.1
  
  Varphi <- exp(varphi0 + t(varphi1)%*%X)
  
  mylist<-list("varphi0"=varphi0,"varphi1"=varphi1,"P.t"=Varphi)
  
  return(mylist)
}


varphi.hat.fast<-function(model_sol,omega,H,x,a,b,X=model_sol$X,t=0){
  dx<-matrix(x-c(0,x[1:length(x)-1]),length(x),1)
  U <- matrix(omega,model_sol$n.X,length(x))
  U <- U + 1i*matrix(a,ncol=1) %*% matrix(x,nrow=1)
  s1<-varphi_multi_omega(model_sol,U,H,X=X,t=t)$P.t
  fx<-outer(x,b,function(r,c)Im(s1[,1]*exp(-1i*r*c))/r)*
    dx[,1]
  varphi.hat<-varphi(model_sol,omega,H,X,t)[[3]][H]/2-1/pi*apply(fx,2,sum)
  return(varphi.hat)
}

varphi.bar.fast<-function(model_sol,omega,H,x,a,b,X=model_sol$X,t=0){
  eps       <-10^-5
  varphi.bar<-(varphi.hat.fast(model_sol,eps*omega,H,x,a,b,X,t)[[1]]-
                 varphi.hat.fast(model_sol,0,H,x,a,b,X,t)[[1]])/eps
  return(varphi.bar)
}


update.model_sol.4.mu_altern <- function(model_sol,mu_altern,
                                         indic_cum = 1,
                                         indic_Euler = 0,
                                         H_if_Euler = model_sol$Tmax
){
  # If indic_cum == 1, then mu_altern defines the growth rate of the
  #     the new variable.
  # If indic_Euler == 1, then the new variable is a price, and 
  #     mu_altern$muprice_0 is updated so that the Euler equation is satisfied.
  # In that case (indic_Euler == 1), mu_altern$muprice_0 are updated until
  #     maturity H_if_Euler
  
  new_model_sol <- model_sol
  n.Z <- model_sol$n.Z
  n.W <- model_sol$n.W
  
  mu_altern_Z <- matrix(mu_altern$muprice_1[1:n.Z],ncol=1)
  mu_altern_W <- matrix(mu_altern$muprice_1[(n.Z+1):(n.Z+n.W)],ncol=1)
  
  indic_altern <- which(model_sol$names.var.X=="Cum_dc")
  
  if(indic_Euler==1){
    # In that case, we'll recompute appropriate values of
    # mu_altern$muprice_0 that satisfy the Euler equations
    mu_altern$muprice_0 <- 0
  }
  
  for(i in 1:(model_sol$Tmax-1)){
    new_model_sol$omega[[i]][indic_altern,] <-
      t(mu_altern_Z) %*% new_model_sol$omega[[i]] + t(mu_altern_W)
    new_model_sol$A1[[i]][indic_altern,]    <- 
      t(mu_altern_Z) %*% new_model_sol$A1[[i]]
    new_model_sol$omega0[[i]][indic_altern] <- 
      ifelse(length(mu_altern$muprice_0)>1,
             t(mu_altern_Z) %*% new_model_sol$omega0[[i]] +
               mu_altern$muprice_0[i],
             t(mu_altern_Z) %*% new_model_sol$omega0[[i]] +
               mu_altern$muprice_0)
    if(indic_cum==1){
      new_model_sol$A1[[i]][indic_altern,indic_altern] <- 1
    }
  }
  new_model_sol$omega.inf[indic_altern,] <-
    t(mu_altern_Z) %*% new_model_sol$omega.inf + t(mu_altern_W)
  new_model_sol$A1.inf[indic_altern,]    <- 
    t(mu_altern_Z) %*% new_model_sol$A1.inf
  new_model_sol$omega0.inf[indic_altern] <- 
    ifelse(length(mu_altern$muprice_0)>1,
           t(mu_altern_Z) %*% new_model_sol$omega0.inf +
             mu_altern$muprice_0[model_sol$Tmax],
           t(mu_altern_Z) %*% new_model_sol$omega0.inf +
             mu_altern$muprice_0)
  if(indic_cum==1){
    new_model_sol$A1.inf[indic_altern,indic_altern] <- 1
  }
  
  
  if(indic_Euler==1){
    # Let's recompute appropriate values of
    # mu_altern$muprice_0 (that satisfy the Euler equations)
    omega <- matrix(0,model_sol$n.X,1)
    omega[indic_altern] <- 1
    prices.ZCRF.bonds   <- varphi(new_model_sol,
                                  omega.varphi = 0*omega,
                                  H = H_if_Euler)
    # Determine sequence of mu_A0t:
    prices.Atilde <- varphi(new_model_sol,
                            omega.varphi = omega,
                            H = H_if_Euler)
    varphi.tilde.A   <-  prices.Atilde$P.t
    varphi.tilde.A_1 <- c(1,varphi.tilde.A[1:(H_if_Euler-1)])
    mu_A0t <- - log(varphi.tilde.A) + log(varphi.tilde.A_1)
    # Update model accordingly:
    mu_A <- mu_altern
    mu_A$muprice_0 <- c(mu_A0t,rep(mu_A0t[H],model$Tmax-H_if_Euler))
    
    new_model_sol <- update.model_sol.4.mu_altern(model_sol,
                                                  mu_altern = mu_A,
                                                  indic_Euler=0)
  }
  
  return(new_model_sol)
}



