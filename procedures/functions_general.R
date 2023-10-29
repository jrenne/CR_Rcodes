#Generalized Functions file
#--------------------------ctrl F: CHANGE

#* CALIBRATION ----
solveParam4c <- function(model){
  model$parameters$A_bar <-
    exp(model$target_vector[["mu_c0"]])/model$parameters$delta -
    1 + model$parameters$delta_K
  model$parameters$sigma_a <- model$target_vector[["sigma_c0"]] *
    (model$parameters$A_bar + 1 - model$parameters$delta_K)
  return(model)
}

solveParam4D <-function(model){
  T.2100  <- c(2,4)
  t.star  <- model$horiz.2100
  
  cst     <-log(model$target_vector["stdCumD4"]^2+model$target_vector["ECumD4"]^2)/
    log(model$target_vector["ECumD4"])
  
  mu_D    <-(2-cst)/(2*cst-2)
  
  ell.D<-matrix(c(0,0),ncol=1)
  
  ell.D<-solve(matrix(c(t.star,t.star,
                        .5*((t.star-1)*T.2100[1]+(t.star+1)*model$vector.ini$ini_Tat),
                        .5*((t.star-1)*T.2100[2]+(t.star+1)*model$vector.ini$ini_Tat)),
                      2,2))%*%
    matrix((-1-mu_D)/(mu_D)*
             c(log(model$target_vector["ECumD2"]),
               log(model$target_vector["ECumD4"])),
           ncol=1)
  
  model$parameters$mu_D <- mu_D
  model$parameters$a_D  <- ell.D[1] * mu_D
  model$parameters$b_D  <- ell.D[2] * mu_D
  
  return(model)
}

solveParam4H <-function(model){
  T.2100 <- c(2,4)
  t.star <- model$horiz.2100
  
  C<-matrix(c(model$target_vector["EH2"]-model$vector.ini$ini_H,
              model$target_vector["EH4"]-model$vector.ini$ini_H),
            ncol=1)
  
  S<-matrix(c(.5*((t.star-1)*T.2100[1]+(t.star+1)*model$vector.ini$ini_Tat),
              .5*((t.star-1)*T.2100[2]+(t.star+1)*model$vector.ini$ini_Tat),
              t.star,t.star),
            2,2)
  SAT <- solve(S)%*%C
  
  model$parameters$b_H <- SAT[1]
  model$parameters$a_H <- SAT[2]
  
  model$parameters$mu_H  <- model$target_vector["stdH4"]^2/
    (2*model$target_vector["EH4"])
  
  return(model)
}

solveParam4N <-function(model,
                        T.2100 = c(2,4),
                        T.2500 = 4){
  
  t.star  <- model$horiz.2100
  
  mu_N  <- model$target_vector["stdCumN4"]^2/(2*model$target_vector["ECumN4"])
  
  kapN  <- matrix(seq(0.6,0.9,length=5000),ncol=1)
  
  C     <- matrix(c(model$target_vector["ECumN2"],
                    model$target_vector["ECumN4"]),
                  ncol=1)
  
  ellN.condkapN<-apply(kapN,1,function(kN){
    N<-matrix(c(rep((1-kN^t.star)/(1-kN),2),
                (1-kN^t.star)/(1-kN)*model$vector.ini$ini_Tat+
                  (T.2100[1]-model$vector.ini$ini_Tat)/t.star*
                  (kN+((t.star-1)*kN-t.star)*kN^t.star)/(1-kN)^2,
                (1-kN^t.star)/(1-kN)*model$vector.ini$ini_Tat+
                  (T.2100[2]-model$vector.ini$ini_Tat)/t.star*
                  (kN+((t.star-1)*kN-t.star)*kN^t.star)/(1-kN)^2),
              2,2)*mu_N
    ellN  <- matrix(c(0,0),ncol=1)
    ellN  <- solve(N)%*%C
    
    return(ellN)})
  
  CumN.inf    <- mu_N/(1-kapN[1:length(kapN)])*
    (t(ellN.condkapN[,1:length(kapN)])%*%
       matrix(c(1,T.2500),ncol=1))
  
  param.CumN  <-rbind(t(kapN),ellN.condkapN,t(CumN.inf),
                      rep(model$target_vector["ECumNinf"],length(kapN)),
                      t(CumN.inf-model$target_vector["ECumNinf"]))
  rownames(param.CumN)<-c("kappaN","ell0.N","ell1.N","CumN.inf",
                          "Target",
                          "Diff")
  
  ellN.kapN  <-param.CumN[,which(abs(param.CumN["Diff",])==
                                   min(abs(param.CumN["Diff",])))]
  
  model$parameters$mu_N    <- mu_N
  model$parameters$kappa_N <- ellN.kapN["kappaN"]
  model$parameters$a_N     <- ellN.kapN["ell0.N"] * mu_N
  model$parameters$b_N     <- ellN.kapN["ell1.N"] * mu_N
  
  return(model)
}


solveParam4T <-function(model){
  # determine mu_T so as to match stdTat2100
  
  print("** Calibration of mu_T by bisection method **")
  
  mu_T.low  <- .000001
  mu_T.high <- .3
  
  model.low  <- model
  model.high <- model
  model.low$parameters$mu_T  <- mu_T.low
  model.high$parameters$mu_T <- mu_T.high
  model_sol.low  <- model_solve(model.low,model$theta0)
  model_sol.high <- model_solve(model.high,model$theta0)
  
  EV_low  <- EV.fct(model_sol.low,model$horiz.2100)
  EV_high <- EV.fct(model_sol.high,model$horiz.2100)
  
  std.T_at.low  <- sqrt(EV_low$VX$T_at[model$horiz.2100])
  std.T_at.high <- sqrt(EV_high$VX$T_at[model$horiz.2100])
  
  # To be sure the higher bound is high enough:
  while(std.T_at.high<model$target_vector["stdTat2100"]){
    model.high <- model
    model.high$parameters$mu_T <- model.high$parameters$mu_T + .2
    model_sol.high <- model_solve(model.high,model$theta0)
    EV_high <- EV.fct(model_sol.high,model$horiz.2100)
    std.T_at.high <- sqrt(EV_high$VX$T_at[model$horiz.2100])
  }
  
  DIFF     <- abs(model$target_vector["stdTat2100"] - std.T_at.low)
  rel.diff <- DIFF/(std.T_at.high - std.T_at.low)
  
  tol  <- .001
  model.interm  <- model
  while(DIFF > tol){
    
    print(paste(" deviation between higher and lower bound: ",toString(std.T_at.high-std.T_at.low),sep=""))
    
    #mu_T.interm <- mu_T.low * (1-rel.diff) + mu_T.high * rel.diff
    mu_T.interm <- mu_T.low * .5 + mu_T.high * .5
    
    model.interm$parameters$mu_T  <- mu_T.interm
    model_sol.interm  <- model_solve(model.interm,model$theta0)
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
  model_new <- model_solve(model,model$theta0)
  
  EV_new <- EV.fct(model_new,model$horiz.2100)
  Var.T_at <- EV_new$VX$T_at[model$horiz.2100]
  
  return(sqrt(Var.T_at)-model$target_vector["stdTat2100"])
}


#* MISCELLANEOUS ----

##----------------Bisection Algorithm
#*bounds is a size 2 vector containing lower and upper bounds respectively
#*bounds[1]<bounds[2]
#*fun is the function to maximize
#*return a vector containing c and fun(c) respectively

bisection <- function(bounds,fun,model){
  a    <- min(bounds)
  b    <- max(bounds)
  tol  <- .001
  i    <- 0
  diff <- 1
  while(diff > tol){
    c   <- (a+b)/2
    f.c <- fun(c,model)
    f.a <- fun(a,model)
    if(f.c*f.a > 0){
      a <- c
    }else{
      b <- c
    }
    i <- i + 1
    diff <- abs(b-a)
    print(diff)
  }
  return(c(c,fun(c)))
}



#parallel calculations for f(a) and f(c)

#----------------Plot with confidence area
make_confidence_area <- function(x,y,pdf_xy,p,tol=10^(-2)){
  # Computes the coordinates of a polygon that delineates a confidence
  #     area associated with probability p.
  # The components of x (and of y) must be equally spaced.
  
  h.x <- x[2] - x[1]
  h.y <- y[2] - y[1]
  
  pdf_xy_h2 <- pdf_xy * h.x * h.y # The sum of all entries of pdf_xy_h2 should ne close to 1.
  
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
                    tstar){
  param <- model$parameters
  a <- tstar * u*param$mu_D/(1-u*param$mu_D) * param$a_D / param$mu_D
  b <- u*param$mu_D/(1-u*param$mu_D) * param$b_D / param$mu_D /2 * ((tstar-1)*Tstar + (tstar+1)*T0)
  return(exp(a+b))
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
                   tstar){
  param <- model$parameters
  a <- u*model$vector.ini$ini_H +
    tstar * u*param$mu_H/(1-u*param$mu_H) * param$a_H/param$mu_H
  b <- u*param$mu_H/(1-u*param$mu_H) * param$b_H/param$mu_H/2 * ((tstar-1)*Tstar + (tstar+1)*T0)
  return(exp(a+b))
}

Param.Gamma0.CumN <- function(model,T0=model$vector.ini$ini_Tat,Tstar,tstar){
  param  <- model$parameters
  kN     <- param$kappa_N
  ell0.N <- param$a_N / param$mu_N
  ell1.N <- param$b_N / param$mu_N
  mu_N   <- param$mu_N
  lambda <- (1-kN^tstar)/(1-kN)*(ell0.N + ell1.N * T0) + 
    ell1.N * (Tstar - T0)/tstar * (kN+((tstar-1)*kN-tstar)*kN^tstar)/(1-kN)^2
  return(list(lambda=lambda,
              mu=mu_N))
}

#------- Fourier transform for permafrost-related emissions:
Fourier.psi <- function(model,gamma,x,
                        T0=model$vector.ini$ini_Tat,
                        Tstar,
                        tstar,psi.LT){
  
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
                             T0=NaN,
                             Tstar,
                             tstar=NaN){
  param <- model$parameters
  aux <- 1/(1-param$kappa_N) * u*param$mu_N/(1-u*param$mu_N) * 
    (param$a_N + param$b_N * Tstar) / param$mu_N
  return(exp(aux))
}


#* SCC ----

#----------SCC, for h>=0
#*if all = TRUE, vector of H \in(0, h).

scc.fct<-function(model_sol,h,
                  all=F,C_0=model_sol$parameters$c0,X=model_sol$X,t=0){  
  mat <- which(model_sol$names.var.X=="M_at")
  mu.u1.c    <-abs(extract(model_sol$mu_u1.t1,mat))
  if(h>(model_sol$Tmax-1)){
    mu.u1.c    <-c(mu.u1.c,rep(abs(model_sol$inf_matx$mu_u1[mat]),
                               h-model_sol$Tmax+1)) 
  }
  omega_c    <-matrix(0,model_sol$n.X,1)
  omega_c[which(model_sol$names.var.X=="Cum_dc")]<-1
  
  if(h==0){
    scc<--C_0*mu_u.t.fct(model_sol)[[2]][which(model_sol$names.var.X=="M_at")]*10^(3)
  } else{
    scc<-multi.lt.fct(model_sol,omega_c,h,X,t)$uX_t.h*
      mu.u1.c[h+1]*C_0*10^(3)
  }
  
  if(all){
    scc.all<--C_0*mu_u.t.fct(model_sol)[[2]][mat]*10^(3)
    for(i in 1:(h-1)){
      scc.all[i+1]<-multi.lt.fct(model_sol,omega_c,i,X)$uX_t.h*
        mu.u1.c[i+1]*C_0*10^(3)
    }
    return(scc.all)/(1-model_sol$parameters$delta)
  }
  return(scc/(1-model_sol$parameters$delta))
}


#* PRICING ----

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
varphi<-function(model_sol,omega.varphi,H,X=model_sol$X,t=0){
  if((t+H)>=(model_sol$Tmax-1)){
    P.pi    <-model_sol$pi
    inf     <-rep(list(model_sol$inf_matx$pi.inf),t+H-(model_sol$Tmax-1)+1)
    P.pi    <-c(P.pi,inf)
    
    P.eta_1 <-model_sol$eta1
    inf     <-rep(list(model_sol$inf_matx$eta1.inf),t+H-(model_sol$Tmax-1)+1)
    P.eta_1 <-c(P.eta_1,inf)
    
    P.eta_0 <-model_sol$eta0
    inf     <-matrix(model_sol$inf_matx$eta0.inf,t+H-(model_sol$Tmax-1)+1,1)
    P.eta_0 <-rbind(P.eta_0,inf)
  }else{
    P.pi    <-model_sol$pi
    P.eta_1 <-model_sol$eta1
    P.eta_0 <-model_sol$eta0
  }
  
  #List of all our U.sh
  U.tk<-list(P.pi[[t+1]]+omega.varphi)
  P.a.pi<-matrix(NaN,H,1)
  for(h in 1:H){
    if((t+h)<=length(model_sol$pi)){
      P.a.pi[h]<-a1.fct(model_sol,P.pi[[t+h]],t+(h-1))
    }else{
      P.a.pi[h]<-a1.fct.inf(model_sol,P.pi[[t+h]]) 
    }
  }
  
  if(H>1){
    for (h in 2:H){
      uk<-matrix(NaN,nrow=model_sol$n.X,ncol=h)
      for(k in 1:(h-1)){
        if((t+k)<(model_sol$Tmax-1)){
          uk[,k] <--P.eta_1[[t+k+1]]-
            b1.fct(model_sol,P.pi[[t+k+1]],t+k)+
            P.pi[[t+k]]
        }else{
          uk[,k] <--P.eta_1[[t+k+1]]-
            b1.fct.inf(model_sol,P.pi[[t+k+1]])+
            P.pi[[t+k]]
        }
      }
      uk[,h]   <- P.pi[[t+h]]+omega.varphi
      U.tk[[h]]<-uk
    }
  }
  
  P.psi<-lapply(1:H,function(h){
    if(h==1){
      multi.lt.fct(model_sol,U.tk[[h]],h,X,t)
    }else{multi.lt.fct.Uh(model_sol,U.tk[[h]],X,t)}
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


#* FOURIER TRANSFORM ----

#------------------------Fourier defined for cdf (u=0)
#*i is the variable in X we are interested in
#*h is the maturity/horizon
#*u is set to 0 to compute the cdf
#*x is a grid for the integral
#*gamma is the set of possible values for the variable i
# psi<-function(model_sol,u,h,i,X=model_sol$X){
#   U    <-matrix(0,model_sol$n.X,length(u))
#   U[i,]<-u
#   psi  <-multi.lt.fct.N(model_sol,U,h,X)
#   return(psi)
# }
# 
# fourier<-function(model_sol,x,gamma,h,i,X=model_sol$X,u=0){
#   dx<-matrix(x-c(0,x[1:length(x)-1]),length(x),1)
#   s1<-psi(model_sol,u+1i*x,h,i,X)
#   fx<-outer(x,gamma,function(r,c)Im(s1[,1]*exp(-1i*r*c))/r)*
#       dx[,1]
#   f <-1/2-1/pi*apply(fx,2,sum)
#   return(f)
# }
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


#* MULTI LT, SIMPLE + COMPLETE ----

#----------------------------------Proposition: Multi-horizon LapT of X
#*U is a matrix of dimension "dim(X)*H", with H>1, the maturity
#*return psi0, psi1 and psi_t.H
multi.lt.fct.Uh<-function(model_sol,Uh,X=model_sol$X,t=0){
  H  <-dim(Uh)[2]
  U  <-list(Uh[,H,drop=F])
  if((t+H)>(model_sol$Tmax-1)){
    a.h<-matrix(c(a1.fct.inf(model_sol,U[[1]]),rep(NaN,H-1)),H,1) 
    i<-1
    #as long as t+H-i is larger than 98=max(a1/b1), we stay time-indep
    while((t+H-i)>(model_sol$Tmax-2)&i<H){
      U[[i+1]]<-Uh[,H-i,drop=F]+b1.fct.inf(model_sol,U[[i]])
      a.h[i+1]<-a1.fct.inf(model_sol,U[[i+1]])
      i<-i+1
    }
    if(i<=(H-1)){
      for (k in i:(H-1)){
        U[[k+1]]<-Uh[,H-k,drop=F]+b1.fct(model_sol,U[[k]],t+H-k)
        a.h[k+1]<-a1.fct(model_sol,U[[k+1]],t+H-(k+1))
      } 
    }
  }else{
    a.h<-matrix(c(a1.fct(model_sol,U[[1]],t+H-1),rep(NaN,H-1)),H,1)   
    for (k in 1:(H-1)){
      U[[k+1]]<-Uh[,H-k,drop=F]+b1.fct(model_sol,U[[k]],t+H-k)
      a.h[k+1]<-a1.fct(model_sol,U[[k+1]],t+H-(k+1))
    }
  }
  if(t>(model_sol$Tmax-2)){
    uX_t.h<-exp(apply(a.h,2,sum)+t(b1.fct.inf(model_sol,U[[H]]))%*%X)
    psi.0 <-apply(a.h,2,sum)
    psi.1 <-b1.fct.inf(model_sol,U[[H]])
  }else{
    uX_t.h<-exp(apply(a.h,2,sum)+t(b1.fct(model_sol,U[[H]],t))%*%X)
    psi.0 <-apply(a.h,2,sum)
    psi.1 <-b1.fct(model_sol,U[[H]],t) 
  }
  mylist<-list("psi.0"=psi.0,"psi.1"=psi.1,"uX_t.h"=uX_t.h)
  return(mylist)
}

#----------------------------------Corollary: Simple Multihorizon LapT
#*h corresponds to the maturity, max 99, then time-independent functions
#*U is a column vector of dimension "dim(X)*1"
#*return psi0, psi1 and psi_t.h
multi.lt.fct<-function(model_sol,U,h,X=model_sol$X,t=0){
  param<-model_sol$parameters
  
  U.h  <-list(U)
  a.sum<-0
  #recursive method with h going to time-indep maturity
  if((t+h)>(model_sol$Tmax-1)){
    i<-1
    while((t+h+1-i)>(model_sol$Tmax-1)){
      a.sum     <-a.sum+a1.fct.inf(model_sol,U.h[[i]])
      U.h[[i+1]]<-b1.fct.inf(model_sol,U.h[[i]])
      
      i<-i+1
    }
    if(i<=h){
      for (k in i:h){
        a.sum     <-a.sum+a1.fct(model_sol,U.h[[k]],t+h-k)
        U.h[[k+1]]<-b1.fct(model_sol,U.h[[k]],t+h-k)
      }  
    }
    uX_t.h<-exp(a.sum+t(U.h[[h+1]])%*%X)
  }else{#Or max of time-dep maturity
    for (k in 1:h){
      a.sum     <-a.sum+a1.fct(model_sol,U.h[[k]],t+h-k)
      U.h[[k+1]]<-b1.fct(model_sol,U.h[[k]],t+h-k)
    }
    uX_t.h<-exp(a.sum+t(U.h[[h+1]])%*%X)
  }
  psi0<-a.sum
  psi1<-U.h[[h+1]]
  
  mylist<-list("psi0"=psi0,"psi1"=psi1,"uX_t.h"=uX_t.h)
  return(mylist)
}

#----------------------------------Corollary .N: Simple Multihorizon LapT
#*U is a matrix of dimension "dim(X)*N"
#*return psi_t.h
# multi.lt.fct.N<-function(model_sol,U,h,X=model_sol$X,t=0){
#   param    <-model_sol$parameters
#   
#   U.h  <-list(U)
#   a.sum<-matrix(0,nrow=dim(U)[2])
#   if((t+h)>(model_sol$Tmax-1)){
#     i<-1
#     while((t+h+1-i)>(model_sol$Tmax-1)){
#       a.sum     <-a.sum+a1.fct.inf(model_sol,U.h[[i]])
#       U.h[[i+1]]<-b1.fct.inf(model_sol,U.h[[i]])
#       
#       i<-i+1
#     }
#     if(i<=h){
#       for (k in i:h){
#         a.sum     <-a.sum+a1.fct(model_sol,U.h[[k]],t+h-k)
#         U.h[[k+1]]<-b1.fct(model_sol,U.h[[k]],t+h-k)
#       } 
#     }
#     uX_t.h<-exp(a.sum+t(U.h[[h+1]])%*%X)
#     
#   }else{
#     for (k in 1:h){
#       a.sum     <-a.sum+a1.fct(model_sol,U.h[[k]],t+h-k)
#       U.h[[k+1]]<-b1.fct(model_sol,U.h[[k]],t+h-k)
#     }
#     uX_t.h<-exp(a.sum+t(U.h[[h+1]])%*%X)
#   }
#   return(uX_t.h)
# }
multi.lt.fct.N<-function(model_sol,U,h,X=model_sol$X,t=0,indic.cum=0){
  param    <-model_sol$parameters
  
  U.h  <-list(U)
  a.sum<-matrix(0,nrow=dim(U)[2])
  if((t+h)>(model_sol$Tmax-1)){
    i<-1
    while((t+h+1-i)>(model_sol$Tmax-1)){
      a.sum     <-a.sum+a1.fct.inf(model_sol,U.h[[i]]+indic.cum*(i>1)*U)
      U.h[[i+1]]<-b1.fct.inf(model_sol,U.h[[i]]+indic.cum*(i>1)*U)
      
      i<-i+1
    }
    if(i<=h){
      for (k in i:h){
        a.sum     <-a.sum+a1.fct(model_sol,U.h[[k]]+indic.cum*(i>1)*U,t+h-k)
        U.h[[k+1]]<-b1.fct(model_sol,U.h[[k]]+indic.cum*(i>1)*U,t+h-k)
      } 
    }
    uX_t.h<-exp(a.sum+t(U.h[[h+1]])%*%X)
    
  }else{
    for (k in 1:h){
      a.sum     <-a.sum+a1.fct(model_sol,U.h[[k]]+indic.cum*(k>1)*U,t+h-k)
      
      U.h[[k+1]]<-b1.fct(model_sol,U.h[[k]]+indic.cum*(k>1)*U,t+h-k)
    }
    uX_t.h<-exp(a.sum+t(U.h[[h+1]])%*%X)
  }
  return(uX_t.h)
}


#* SOLVING SDF ----

#------------------Proposition: Function that returns mu_u1 and mu_u0 @0 
#*Tend: General case with Tmax=100, end up with t=2020
#*If want a specific date--> mu_u.t.fct.all
#*return mu_u0.t, mu_u1.t
mu_u.t.fct<-function(model_sol,Tend=model_sol$Tmax){
  param    <-model_sol$parameters
  mu_u1.inf<-model_sol$mu_u1
  mu_u0.inf<-model_sol$mu_u0
  
  mu_c1  <-model_sol$mu_c1
  mu_c0  <-model_sol$mu_c0
  mu_u1.t<-mu_u1.inf                                                            #mu_u1.99
  mu_u0.t<-mu_u0.inf                                                            #mu_u0.99
  #mu_u0 and mu_u1.98-0,99=inf.
  for (i in 2:Tend){
    b.t    <-b1.fct(model_sol,(1-param$gamma)*(mu_u1.t+mu_c1),model_sol$Tmax-i)
    a.t    <-a1.fct(model_sol,(1-param$gamma)*(mu_u1.t+mu_c1),model_sol$Tmax-i)
    mu_u1.t<-param$delta/(1-param$gamma)*b.t
    mu_u0.t<-param$delta*(mu_u0.t+mu_c0)+
      param$delta/(1-param$gamma)*a.t
  }
  
  mylist<-list("mu_u0.t"=mu_u0.t,"mu_u1.t"=mu_u1.t)
  return(mylist)
}

#------------------list with all mu_u1.t and mu_u0.t
#*from t=1 to t=98, inf(t=99) not in list
mu_u.t.fct.all<-function(model_sol){
  param    <-model_sol$parameters
  mu_u1.inf<-model_sol$mu_u1
  mu_u0.inf<-model_sol$mu_u0
  Tmax     <-model_sol$Tmax
  
  mu_c1  <-model_sol$mu_c1
  mu_c0  <-model_sol$mu_c0
  
  mu_u1.t        <-list()
  mu_u1.t[[Tmax]]<-mu_u1.inf
  mu_u0.t        <-matrix(NaN,Tmax,1)
  mu_u0.t[Tmax]  <-mu_u0.inf
  for (i in 1:(Tmax-1)){
    mu_u1.t[[Tmax-i]]<-param$delta/(1-param$gamma)*
      b1.fct(model_sol,
             (1-param$gamma)*(mu_u1.t[[Tmax-i+1]]+mu_c1),
             Tmax-i-1)
    mu_u0.t[Tmax-i]  <-param$delta*(mu_u0.t[Tmax-i+1]+mu_c0)+
      param$delta/(1-param$gamma)*
      a1.fct(model_sol,
             (1-param$gamma)*(mu_u1.t[[Tmax-i+1]]+mu_c1),
             Tmax-i-1)
  }
  #no inf
  mu_u1.t<-mu_u1.t[-Tmax]
  mu_u0.t<-mu_u0.t[-Tmax]
  #no 2020
  mu_u1  <-mu_u1.t[-1]
  mu_u0  <-mu_u0.t[-1]
  vec_date<-seq(model_sol$vec_date[2],by=model_sol$tstep,length=Tmax-2)         #infinity not counted
  
  mylist<-list("mu_u0.t1"=mu_u0,"mu_u1.t1"=mu_u1,"date"=vec_date)
  return(mylist)
}


#* SOLVE MODEL, MEAN AND SIMUL ----

#----------------------------------Solve model for infinite mitig mu=1
#*theta is a set of ini cond. for the optim. of mitig rate mu
model_solve<-function(model,theta,
                      indic_mitig=TRUE,mu.chosen=rep(param$mu0,model$Tmax),
                      mu_altern=model$Cum_dc){
  #Data preparation
  model_sol<-model
  param    <-model$parameters
  tstep    <-model$tstep
  GN       <-0
  Tmax     <-model$Tmax
  
  # Create indicators of position of variables (Variables in Z):
  for(i in 1:model$n.Z){
    eval(parse(text = gsub(" "," ",
                           paste("indic.",model$names.var.X[i],
                                 "<- which(model$names.var.X=='",model$names.var.X[i],
                                 "')",sep=""))))}
  
  # Create indicators of position of variables (Variables in W):
  for(i in (model$n.Z+1):(model$n.Z+model$n.W)){
    eval(parse(text = gsub(" "," ",
                           paste("indic.",model$names.var.X[i],
                                 "<- which(model$names.var.X=='",model$names.var.X[i],
                                 "')-model$n.Z",sep=""))))}
  
  ell1.D   <-matrix(0,nrow=model_sol$n.Z+model_sol$n.W)
  ell1.D[indic.T_at]   <- model_sol$parameters$b_D/model_sol$parameters$mu_D
  model_sol[["ell1.D"]]<- ell1.D
  
  ell1.N   <-matrix(0,nrow=model_sol$n.Z+model_sol$n.W)
  ell1.N[indic.T_at]   <- model_sol$parameters$b_N/model_sol$parameters$mu_N
  model_sol[["ell1.N"]]<- ell1.N
  
  ell1.T   <-matrix(0,nrow=model_sol$n.Z+model_sol$n.W)
  ell1.T[indic.Forc]   <- param$xi_1
  ell1.T[indic.T_at]   <- 1-param$xi_1*(param$tau/param$nu + param$xi_2)
  ell1.T[indic.T_lo]   <- param$xi_1*param$xi_2
  ell1.T               <- ell1.T/param$mu_T
  model_sol[["ell1.T"]]<- ell1.T
  
  ell1.H   <-matrix(0,nrow=model_sol$n.Z+model_sol$n.W)
  ell1.H[indic.T_at]    <- param$b_H/param$mu_H
  model_sol[["ell1.H"]] <- ell1.H
  
  #Determine the matrix of shocks                                               #+CHANGE in mu_dep
  omega.star.inf              <- matrix(0,model_sol$n.Z,model_sol$n.W)
  omega.star.inf[c(indic.delc,indic.y_tilde),indic.eta_A] <-
    param$sigma_a/(param$A_bar+1-param$delta_K)
  #omega.star.inf[indic.E_ind,indic.eta_E] <- 0*param$sigma_E
  #omega.star.inf[indic.Forc,indic.eta_F]  <- param$sigma_F
  omega.star.inf[indic.Cum_dc,] <- t(mu_altern$muprice_1[(model_sol$n.Z+1):
                                                           (model_sol$n.Z+model_sol$n.W)])
  
  omega.star.inf[indic.delc,indic.D]     <- -1 #shock D
  omega.star.inf[indic.Cum_D,indic.D]    <- -1
  omega.star.inf[indic.E,indic.N]        <-  1 #shock N
  omega.star.inf[indic.T_at,indic.T_atW] <-  1 
  omega.star.inf[indic.H,indic.HW]       <-  1 
  
  model_sol[["omega.star.inf"]]<-omega.star.inf
  
  ##############################################################################
  #Long term matrices and Gauss Newton Algorithm
  ##############################################################################
  #CHANGE vectorial rep.inf.
  A0.star.inf       <- diag(model_sol$n.Z)
  A0.star.inf[indic.delc,indic.H]    <- param$b_sk
  A0.star.inf[indic.E,indic.E_ind]   <- -1
  A0.star.inf[indic.Forc,indic.M_at] <- -param$tau/(log(2)*param$m_pi*param$m0)
  A0.star.inf[indic.Cum_E,indic.E]   <- -1
  A0.star.inf[indic.Cum_dc,]         <- -t(mu_altern$muprice_1[1:model_sol$n.Z])
  A0.star.inf[indic.Cum_dc,indic.Cum_dc] <- 1
  
  varphi           <- matrix(0,3,3)
  varphi[1,1]      <- param$varphi_11
  varphi[2,1]      <- param$varphi_12
  varphi[1,2]      <- param$varphi_21
  varphi[2,2]      <- param$varphi_22
  varphi[3,2]      <- param$varphi_23
  varphi[2,3]      <- param$varphi_32
  varphi[3,3]      <- param$varphi_33
  
  A1.star.inf           <- matrix(0,model_sol$n.Z,model_sol$n.Z)
  A1.star.inf[indic.delc,indic.H] <- param$b_sk
  A1.star.inf[indic.y_tilde,indic.y_tilde]    <- 1
  A1.star.inf[indic.E_ind,indic.y_tilde]      <- 0
  A1.star.inf[indic.M_at,indic.E] <- model_sol$tstep/3.666
  A1.star.inf[indic.M_at:indic.M_lo,indic.M_at:indic.M_lo]  <- varphi%^%(model_sol$tstep)                           
  A1.star.inf[indic.T_lo,indic.T_at]    <- param$xi_3
  A1.star.inf[indic.T_lo,indic.T_lo]    <- 1-param$xi_3
  A1.star.inf[indic.Cum_D,indic.Cum_D]  <- 1
  A1.star.inf[indic.Cum_E,indic.Cum_E]  <- 1
  A1.star.inf[indic.Cum_dc,indic.Cum_dc]<- 1
  A1.star.inf[indic.H,indic.H]          <- 1
  
  omega0.star.inf      <- matrix(0,model_sol$n.Z,1)
  omega0.star.inf[indic.delc,1]  <- log(param$delta)+log(param$A_bar+1-param$delta_K)
  omega0.star.inf[indic.E,1]     <- 0
  omega0.star.inf[indic.E_ind,1] <- 0
  omega0.star.inf[indic.Forc,1]  <- param$tau/log(2)*(log(param$m0)-1)+param$phi_1
  if(length(mu_altern$muprice_0)==1){
    omega0.star.inf[indic.Cum_dc,1] <- mu_altern$muprice_0
  }else{
    omega0.star.inf[indic.Cum_dc,1] <- mu_altern$muprice_0[model$Tmax]
  }
  
  A1.inf      <-solve(A0.star.inf)%*%A1.star.inf
  omega0.inf  <-solve(A0.star.inf)%*%omega0.star.inf
  omega.inf   <-solve(A0.star.inf)%*%omega.star.inf
  
  model_sol[["varphi"]]       <-varphi
  model_sol[["A0.star.inf"]]  <-A0.star.inf
  model_sol[["A1.inf"]]       <-A1.inf
  model_sol[["omega0.inf"]]   <-omega0.inf
  model_sol[["omega.inf"]]    <-omega.inf
  
  ##Infinite components (t+100)
  # solution initial values
  betawGN<- cbind(matrix(0,model_sol$n.Z+model_sol$n.W,model_sol$n.Z),          #CHANGE IF ADD GAMMA0
                  rbind(
                    matrix(0,model_sol$n.Z,model_sol$n.eta),
                    t(model_sol$parameters$Phi),
                    matrix(0,model_sol$n.W-model_sol$n.eta,model_sol$n.eta)
                  ),
                  model_sol$parameters$mu_D*model_sol$ell1.D,
                  model_sol$parameters$mu_N*model_sol$ell1.N,
                  model_sol$parameters$mu_T*model_sol$ell1.T,
                  model_sol$parameters$mu_H*model_sol$ell1.H
  )
  betaGN <- cbind(rbind(t(model_sol$A1.inf),
                        matrix(0,model_sol$n.W,model_sol$n.Z)),
                  matrix(0,model_sol$n.Z+model_sol$n.W,model_sol$n.W)
  )+
    betawGN %*% rbind(matrix(0,model_sol$n.Z,model_sol$n.Z+model_sol$n.W),
                      cbind(t(model_sol$omega.inf),diag(model_sol$n.W))
    )
  model_sol[["M"]]<-betaGN
  muc1GN  <- rbind(1,matrix(0,model_sol$n.Z+model_sol$n.W-1))
  muu1_ini<- model_sol$parameters$delta*
    solve(diag(model_sol$n.Z+model_sol$n.W)-model_sol$parameters$delta*betaGN)%*%
    betaGN%*%muc1GN
  GN      <-muu1_ini
  
  Newton<-GN.function(model_sol,GN)
  
  mu_c1 <-Newton$mu_c1
  mu_c0 <-0
  mu_u1 <-Newton$mu_u1
  mu_u0 <-param$delta/(1-param$delta)*
    mu_c0+param$delta/((1-param$delta)*(1-param$gamma))*
    a1.fct.inf(model_sol,(1-param$gamma)*(mu_u1+mu_c1))
  
  model_sol[["mu_c1"]]<-mu_c1
  model_sol[["mu_c0"]]<-mu_c0
  model_sol[["mu_u1"]]<-mu_u1
  model_sol[["mu_u0"]]<-mu_u0
  model_sol[["ite"]]  <-Newton$ite
  model_sol[["dev"]]  <-Newton$dev
  
  ##############################################################################
  #Optimization of mu - X and exogenous vectors
  ##############################################################################
  #####
  #Exogenous variables
  #####
  f_ex  <-matrix(NaN, nrow=Tmax, 1)
  E_land<-matrix(NaN, nrow=Tmax, 1)
  gsigma<-matrix(NaN, nrow=Tmax, 1)
  sigma <-matrix(NaN, nrow=Tmax, 1)
  mu    <-matrix(NaN, nrow=Tmax, 1)
  bp    <-matrix(NaN, nrow=Tmax, 1)
  bc    <-matrix(NaN, nrow=Tmax, 1)
  
  #Radiative forcings
  f_ex               <- matrix(rep(param$phi_0,Tmax),Tmax,1)                        
  f_ex[1:17]         <- f_ex[1:17]+(1/17)*(param$phi_1-param$phi_0)*((1:17)-1)
  f_ex[18:Tmax]      <- f_ex[18:Tmax] + (param$phi_1-param$phi_0)
  
  #Emissions from deforestation
  E_land[1:Tmax]<- param$eps_0*(1-param$rho)**((1:(Tmax))-1)           
  
  #Carbon intensity
  gsigma[1]<-param$gsigma1                                                      
  for(i in 2:Tmax) gsigma[i] <- gsigma[i-1]*((1+param$delsigma)**tstep)    
  sigma[1] <- param$sigma0                                                         
  for(i in 2:Tmax) sigma[i]  <- sigma[i-1] * exp(gsigma[i-1] * tstep)      
  
  #Abatement cost exogenous components
  bp[1:Tmax]<-param$pback*(1-param$gback)**((1:Tmax)-1)              
  bc[1:Tmax]<-bp[1:Tmax]*sigma[1:Tmax]/param$theta2/1000
  
  #####
  #Initial states in 2020                                                       #CHANGE STATES INI COND.
  #####
  Z    <- matrix(NaN,model_sol$n.Z,1)
  Z[indic.delc] <- log(param$delta)+log((1-bc[1]*param$mu0**(param$theta2))*param$A_bar+
                                          1-param$delta_K)
  Z[indic.y_tilde] <- model_sol$vector.ini$ini_tildey
  Z[indic.E] <- model_sol$vector.ini$ini_E
  Z[indic.E_ind] <- model_sol$vector.ini$ini_Eind
  Z[indic.Forc]  <- model_sol$vector.ini$ini_F
  Z[indic.M_at]  <- model_sol$vector.ini$ini_Mat
  Z[indic.M_up]  <- model_sol$vector.ini$ini_Mup
  Z[indic.M_lo]  <- model_sol$vector.ini$ini_Mlo
  Z[indic.T_at]  <- model_sol$vector.ini$ini_Tat 
  Z[indic.T_lo]  <- model_sol$vector.ini$ini_Tlo
  Z[indic.Cum_D] <- model_sol$vector.ini$ini_CumD
  Z[indic.Cum_E] <- model_sol$vector.ini$ini_CumE
  Z[indic.Cum_dc]<- model_sol$vector.ini$ini_Cumdelc
  Z[indic.H]     <- model_sol$vector.ini$ini_H
  W    <- matrix(0,model_sol$n.W,1)
  #t=0, value of the state variables at time 0
  #X=[delc,tilde_y,E,E_ind,Forc,M_at,M_up,M_lo,T_at,T_lo,Cum_D,Cum_E,Cum_dc,W]
  X<-rbind(Z,W)
  
  model_sol$vector.ini$ini_delc    <- Z[indic.delc]
  model_sol$vector.ini$ini_Cumdelc <- Z[indic.delc]
  
  #Keep all elements st [1] is t for t=0
  model_sol[["f_ex"]]  <-f_ex
  model_sol[["E_land"]]<-E_land
  model_sol[["gsigma"]]<-gsigma
  model_sol[["sigma"]] <-sigma
  model_sol[["bp"]]    <-bp
  model_sol[["bc"]]    <-bc
  
  model_sol[["X"]]     <-X
  model_sol[["n.X"]]   <-length(X)
  
  #####
  #Optimization calculations
  #####
  if (indic_mitig){
    opt<-list()
    for (i in 1:length(theta)){
      opt[[i]] <- res.optim(model_sol,theta[[i]],Tend=Tmax,X=X)
    }
    maxi     <-min(unlist(extract(opt,2)))
    best     <-which(maxi==unlist(extract(opt,2)))[1]
    theta.opt<-unlist(opt[[best]][1])
    mu       <-mu.function(model_sol,theta.opt)
    
    model_sol[["theta.opt"]]<-theta.opt
    model_sol[["u0"]]       <-abs(maxi)
  }else{
    mu   <- mu.chosen
  }
  
  model_sol[["mu"]]   <-mu
  #Date 0 storage
  list2015             <-list()
  
  ##############################################################################
  #Exogenous variables depending on mu optimized (or chosen)
  #####
  model_matrix<-mu_dep(model_sol,model_sol$mu,mu_altern=mu_altern)
  #####
  
  #Keep all elements st [1] is t for t=0
  model_sol[["AC"]]      <-model_matrix$AC
  model_sol[["lambda"]]  <-model_matrix$lambda
  
  #remove first period st starting from t=0, first element is t=1
  model_sol[["A1"]]      <-model_matrix$A1[-1]
  model_sol[["omega"]]   <-model_matrix$omega[-1]
  model_sol[["omega0"]]  <-model_matrix$omega0[-1]
  
  #Date 0 storage
  list2015[["A1"]]       <-model_matrix$A1[[1]]
  list2015[["omega"]]    <-model_matrix$omega[[1]]
  list2015[["omega0"]]   <-model_matrix$omega0[[1]]
  
  if (indic_mitig==FALSE){
    mu_u1.t  <-model_sol$mu_u1
    mu_u0.t  <-model_sol$mu_u0
    mu_c1    <-model_sol$mu_c1
    mu_c0    <-model_sol$mu_c0
    
    mu_u.1<-mu_u.t.fct(model_sol)
    
    if(is.na(mu_u.1[[1]])){
      u0 <- 10000
      model_sol[["u0"]] <-u0
    }else{
      u0<-log(model_sol$parameters$c0)+mu_u.1[[1]]+t(mu_u.1[[2]])%*%model_sol$X
      model_sol[["u0"]] <-u0
    }
  }
  
  
  #####
  #*construct mu_u1.t
  #####
  #Proposition 9: SDF construction
  P.mu_u.t1 <-mu_u.t.fct.all(model_sol)
  P.mu_u1.t1<-P.mu_u.t1[[2]]
  
  P.pi  <-lapply(1:(Tmax-2),                                                    
                 function(t){(1-param$gamma)*P.mu_u1.t1[[t]]-
                     param$gamma*mu_c1})
  
  #####
  #*construct eta.0                                                             
  #####
  P.eta_0  <-matrix(unlist(lapply(1:(length(P.pi)),function(t){
    -log(param$delta)+
      a1.fct(model_sol,(1-param$gamma)*(P.mu_u1.t1[[t]]+mu_c1),t-1)-
      a1.fct(model_sol,P.pi[[t]],t-1)})),                             #t=0 for a.t
    length(P.pi),1)
  
  #####
  #*construct eta.1                                                             
  #####
  P.eta_1  <-lapply(1:(length(P.pi)),
                    function(t){b1.fct(model_sol,
                                       (1-param$gamma)*(P.mu_u1.t1[[t]]+mu_c1),
                                       t-1)-
                        b1.fct(model_sol,P.pi[[t]],t-1)})               #t=0 for b.t
  
  
  model_sol[["mu_u1.t1"]]<-P.mu_u1.t1
  model_sol[["pi"]]      <-P.pi
  model_sol[["eta1"]]    <-P.eta_1
  model_sol[["eta0"]]    <-P.eta_0
  
  #Infinite economy:
  inf_matx<-list("mu_u1"=model_sol$mu_u1)
  
  inf_matx[["pi.inf"]]  <-(1-param$gamma)*inf_matx$mu_u1-param$gamma*mu_c1
  inf_matx[["eta0.inf"]]<--log(param$delta)+
    a1.fct.inf(model_sol,(1-param$gamma)*(inf_matx$mu_u1+mu_c1))-
    a1.fct.inf(model_sol,inf_matx$pi.inf)
  inf_matx[["eta1.inf"]]<-b1.fct.inf(model_sol,
                                     (1-param$gamma)*(inf_matx$mu_u1+mu_c1))-
    b1.fct.inf(model_sol,inf_matx$pi.inf)
  
  
  
  model_sol$ini_matx<-c(list2015,model_sol$ini_matx)
  model_sol$inf_matx<-c(inf_matx,model_sol$inf_matx)
  ##############################################################################
  return(model_sol)
}


#* Functions mitig mu optim ----
#-------------------------------Prevent agents to mitigate too quickly
mu.function <- function(model_sol,theta,
                        date.min.mu.equal.1=0,max.mu.date.0=.2){
  a <- -log(max.mu.date.0) + abs(theta[1])
  b <- abs(theta[2])
  
  if(a<b*date.min.mu.equal.1){
    b <- a/date.min.mu.equal.1
  }
  mu  <-pmin(exp(-a+b*(1:model_sol$Tmax)),1)
  return(mu)
}


#-------------------------------Compute the vectorial repr. of the model + exogenous dependent fct(mu)
mu_dep<-function(model_sol,mu,mu_altern=model_sol$Cum_dc){
  param     <-model_sol$parameters
  Tmax      <-model_sol$Tmax
  
  AC     <-matrix(NaN, nrow=Tmax, 1)
  lambda <-matrix(NaN, nrow=Tmax, 1)
  
  #Exogenous equations independent of mu
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
  
  #Exogenous equations dependent on mu
  #Abatement Cost
  AC[1]        <-bc[1]*param$mu0**(param$theta2)#bc[1]*mu[1]**(param$theta2)#                                   
  AC[2:Tmax]   <-bc[2:Tmax]*mu[2:Tmax]**param$theta2  
  
  #lambda
  lambda[1]  <- sigma[1]*(1-param$mu0)*param$q0*#sigma[1]*(1-mu[1])*param$q0*#
    exp((log(param$delta)+log((1-AC[1])*param$A_bar+1-param$delta_K)+
           1/2*((1-AC[1])*param$sigma_a/
                  ((1-AC[1])*param$A_bar+1-param$delta_K))^2))
  
  for (i in 2:Tmax){
    lambda[i] <-(1-mu[i])*sigma[i]*param$q0*
      exp(sum(log(param$delta)+log((1-AC[1:i])*param$A_bar+1-param$delta_K)+
                1/2*((1-AC[1:i])*param$sigma_a/
                       ((1-AC[1:i])*param$A_bar+1-param$delta_K))^2))
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
  
  A0.star<-model_sol$A0.star.inf
  
  A1.star<-list()
  for (i in 1:Tmax){
    A1_i <- matrix(0,model_sol$n.Z,model_sol$n.Z)
    A1_i[indic.y_tilde,indic.y_tilde] <- 1
    A1_i[indic.E_ind,indic.y_tilde]<- lambda[i]
    A1_i[indic.M_at,indic.E]       <- model_sol$tstep/3.666
    A1_i[indic.M_at:indic.M_lo,indic.M_at:indic.M_lo]  <- model_sol$varphi%^%(model_sol$tstep)
    A1_i[indic.T_lo,indic.T_at]    <- param$xi_3
    A1_i[indic.T_lo,indic.T_lo]    <- 1-param$xi_3
    A1_i[indic.Cum_D,indic.Cum_D]  <- 1
    A1_i[indic.Cum_E,indic.Cum_E]  <- 1
    A1_i[indic.Cum_dc,indic.Cum_dc]<- 1
    A1_i[indic.H,indic.H]          <- 1
    A1_i[indic.delc,indic.H]       <-  param$b_sk
    A1.star[[i]]   <- A1_i 
  }
  
  omega0.star<-list()
  for(i in 1:Tmax){
    omega0_i        <- matrix(0,model_sol$n.Z,1)
    if(sum(model_sol$mu_c)==0){
      omega0_i[indic.delc,1]   <- log(param$delta)+
        log((1-AC[i])*param$A_bar+1-param$delta_K)
    }else{# in that case, we want to impose a growth trajectory
      omega0_i[indic.delc,1]   <- model_sol$mu_c[i]
    }
    omega0_i[indic.E,1]     <- E_land[i]
    omega0_i[indic.E_ind,1] <- lambda[i]
    omega0_i[indic.Forc,1]  <- param$tau/log(2)*(log(param$m0)-1)+f_ex[i]
    if(length(mu_altern$muprice_0)==1){
      omega0_i[indic.Cum_dc,1]  <- mu_altern$muprice_0
    }else{
      omega0_i[indic.Cum_dc,1]  <- mu_altern$muprice_0[i]
    }
    omega0.star[[i]]<- omega0_i
  }
  
  omega.star<-list()
  for(i in 1:Tmax){
    omega_i <- matrix(0,model_sol$n.Z,model_sol$n.W)
    omega_i[c(indic.delc,indic.y_tilde),indic.eta_A] <- param$sigma_a*
      (1-AC[i])/((1-AC[i])*param$A_bar+1-param$delta_K)
    #omega_i[indic.E_ind,indic.eta_E] <- lambda[i]*param$sigma_E
    #omega_i[indic.Forc,indic.eta_F]  <- param$sigma_F
    omega_i[indic.Cum_dc,] <- t(mu_altern$muprice_1[(model_sol$n.Z+1):
                                                      (model_sol$n.Z+model_sol$n.W)])
    
    omega_i[indic.delc,indic.D]     <- -1 #shock D
    omega_i[indic.Cum_D,indic.D]    <- -1 #shock D
    omega_i[indic.E,indic.N]        <-  1 #shock N on E
    omega_i[indic.T_at,indic.T_atW] <-  1 
    omega_i[indic.H,indic.HW]       <-  1
    
    omega.star[[i]]      <- omega_i
  }
  A1      <-lapply(1:Tmax,function(i)solve(A0.star)%*%
                     A1.star[[i]])
  omega0  <-lapply(1:Tmax,function(i)solve(A0.star)%*%
                     omega0.star[[i]])
  omega   <-lapply(1:Tmax,function(i)solve(A0.star)%*%
                     omega.star[[i]]) 
  # }
  
  mylist  <-list("A1"=A1,"omega0"=omega0,"omega"=omega,
                 "AC"=AC,"lambda"=lambda)
  return(mylist)
}


#-------------------------------Results of optimization problem of mu
#*linked to utility function
#*return the value of optimization and the a/b (theta1/2)

res.optim <-function(model_sol,theta,Tend=model_sol$Tmax,X=model_sol$X){
  res.optim<-optim(par=theta,utility.optim,
                   model_sol=model_sol,
                   Tend=Tend,
                   X=X,
                   gr = NULL,
                   method="Nelder-Mead",
                   #method="CG",
                   #method="BFGS",
                   control=list(trace=F,maxit=model_sol$MAXIT),
                   hessian=FALSE)
  mylist<-list("res.optim$par"=res.optim$par,
               "res.optim$value"=res.optim$value)
  return(mylist)
}

#-------------------------------construction of u0 for optim of mitig rate mu
utility.optim<-function(model_sol,theta,Tend=model_sol$Tmax,X=model_sol$X){
  param     <-model_sol$parameters
  tstep     <-model_sol$tstep
  omega.star<-model_sol$omega.star
  Tmax      <-model_sol$Tmax
  
  #utility
  mu_u1.inf<-model_sol$mu_u1
  mu_u0.inf<-model_sol$mu_u0
  mu_c1    <-model_sol$mu_c1
  mu_c0    <-model_sol$mu_c0                                                   
  mu_u1.t  <-mu_u1.inf                                                            
  mu_u0.t  <-mu_u0.inf
  
  #Parameters to optimize
  a <-theta[1]
  b <-theta[2]
  
  #List of variables
  mu<-matrix(NaN, nrow=Tmax, 1)
  
  #Emissions control rate
  mu[1:Tmax]<-mu.function(model_sol,c(a,b))
  
  model_matrix<-mu_dep(model_sol,mu)
  omega0      <-model_matrix$omega0
  omega       <-model_matrix$omega
  A1          <-model_matrix$A1
  
  model_sol[["A1"]]    <-A1[-1]
  model_sol[["omega"]] <-omega[-1]
  model_sol[["omega0"]]<-omega0[-1]
  
  mu_u.1<-mu_u.t.fct(model_sol,Tend) 
  #if NaN, returns high utility s.t. I don't mitigate and big loss
  if(is.na(mu_u.1[[1]])){
    return(-10000)
  }
  
  u0<-log(param$c0)+mu_u.1[[1]]+t(mu_u.1[[2]])%*%X
  return(-u0)
}



#* Function for theoretical mean and variance ----
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
    omega0<-model_sol$omega0
    inf   <-rep(list(model_sol$omega0.inf),t-model_sol$Tmax+1)
    omega0<-c(omega0,inf)
    
    omega <-model_sol$omega
    inf   <-rep(list(model_sol$omega.inf),t-model_sol$Tmax+1)
    omega <-c(omega,inf)
    
    A1    <-model_sol$A1
    inf   <-rep(list(model_sol$A1.inf),t-model_sol$Tmax+1)
    A1    <-c(A1,inf)
  }else{
    omega0<-model_sol$omega0
    omega <-model_sol$omega
    A1    <-model_sol$A1
  }
  
  X     <-model_sol$X
  n.Z   <-model_sol$n.Z
  n.W   <-model_sol$n.W
  n.eta <-model_sol$n.eta
  
  #Proposition 4+5: Conditional mean and variance of W
  #alpha1.w, beta1.w, alpha2.w, beta2.w                                         #CHANGE if add gamma0
  alpha1.w    <-lapply(1:t,function(i){
    rbind(matrix(0,n.Z+n.eta,1),
          param$a_D,
          if(i<model_sol$Tmax){
            param$a_N*param$kappa_N^i
          } else{
            param$a_N*0
          },
          0, # T_at
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
      0, # T_at
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
  
  beta2.w.T                      <-matrix(0,n.W,n.W)
  beta2.w.T[(n.eta+3),(n.eta+3)] <-2*param$mu_T^2
  
  beta2.w.H                      <-matrix(0,n.W,n.W)
  beta2.w.H[(n.eta+4),(n.eta+4)] <-2*param$mu_H^2
  
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
  
  #Conditional mean for all h
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
  
  EX<-lapply(1:model_sol$n.X,function(i){extract(EXh,i)})                       #CHANGE if add var/shock
  
  names(EX) <- model_sol$names.var.X
  
  EV[["EX"]]<-EX
  
  #Conditional Variance#
  
  
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
  
  
  #VXh#
  
  
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
  
  
  
  #vecVXh
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


#* Simulations ----
#CHANGE
#------------------Function to simulate state variables of our model
#*New approximation of the radiative forcings (linear)
#*nb.simul is number of periods for 1 simulation \in(1,99)
#*nb.traj is the number of simulations done for nb.simul periods
#*for t>99, time-independent matrices
simul.function<-function(model_sol,nb.simul.t,nb.traj){
  #useful parameters
  nb.simul<-nb.simul.t+1 #take into account the t=0
  
  tstep<-model_sol$tstep
  t    <-1:nb.simul
  param<-model_sol$parameters
  
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
  
  #Shocks
  eta  <-matrix(0, nrow=model_sol$n.eta, ncol=nb.traj)
  D    <-matrix(0, nrow=nb.simul, ncol=nb.traj)
  N    <-matrix(0, nrow=nb.simul, ncol=nb.traj)
  T_atW<-matrix(0, nrow=nb.simul, ncol=nb.traj)
  HW   <-matrix(0, nrow=nb.simul, ncol=nb.traj)
  
  # Create indicators of position of variables (Variables in Z):
  for(i in 1:model$n.Z){
    eval(parse(text = gsub(" "," ",
                           paste("indic.",model$names.var.X[i],
                                 "<- which(model$names.var.X=='",model$names.var.X[i],
                                 "')",sep=""))))}
  
  # Create indicators of position of variables (Variables in W):
  for(i in (model$n.Z+1):(model$n.Z+model$n.W)){
    eval(parse(text = gsub(" "," ",
                           paste("indic.",model$names.var.X[i],
                                 "<- which(model$names.var.X=='",model$names.var.X[i],
                                 "')-model$n.Z",sep=""))))}
  
  ##endogenous equations
  #Initial conditions
  delc   <-model_sol$X[indic.delc]                                                       
  y_tilde<-model_sol$X[indic.y_tilde]
  E      <-model_sol$X[indic.E]                                                       #CDICE2021
  E_ind  <-model_sol$X[indic.E_ind]                                                       #CDICE2021
  Forc   <-model_sol$X[indic.Forc]                                                       #CDICE2021
  M_at   <-model_sol$X[indic.M_at]                                                       #CDICE2021,mateq0
  M_up   <-model_sol$X[indic.M_up]                                                       #CDICE2021,mueq0
  M_lo   <-model_sol$X[indic.M_lo]                                                       #CDICE2021,mleq0
  T_at   <-model_sol$X[indic.T_at]                                                       #CDICE2021,tat0
  T_lo   <-model_sol$X[indic.T_lo]                                                      #CDICE2021,tlo0
  Cum_D  <-model_sol$X[indic.Cum_D]
  Cum_E  <-model_sol$X[indic.Cum_E]
  Cum_dc <-model_sol$X[indic.Cum_dc]
  H      <-model_sol$X[indic.H]
  
  
  Z<-list(matrix(
    rep(c(delc,y_tilde,E,E_ind,Forc,M_at,M_up,M_lo,T_at,T_lo,Cum_D,Cum_E,
          Cum_dc,H),
        nb.traj),
    model_sol$n.Z,nb.traj))
  
  Z[[nb.simul+1]]<-D                                                            #CHANGE if gamma0 not only on Tat
  Z[[nb.simul+2]]<-N
  Z[[nb.simul+3]]<-T_atW
  Z[[nb.simul+4]]<-HW
  
  X<-list(model_sol$X)
  
  for (i in 2:nb.simul) {
    #Shocks
    Z[[nb.simul+1]][i,]    <-
      rgamma(nb.traj,rpois(nb.traj,
                           pmax(0,param$a_D/param$mu_D+
                                  t(model_sol$ell1.D[1:model_sol$n.Z,drop=FALSE])%*%
                                  Z[[i-1]])),scale=param$mu_D)
    
    Z[[nb.simul+2]][i,]    <-
      rgamma(nb.traj,rpois(nb.traj,pmax(0,param$kappa_N^(i-1)*
                                          (param$a_N/param$mu_N+
                                             t(model_sol$ell1.N[1:model_sol$n.Z,drop=FALSE])%*%
                                             Z[[i-1]]))),scale=param$mu_N)
    
    Z[[nb.simul+3]][i,]    <-
      rgamma(nb.traj,rpois(nb.traj,pmax(0,0+
                                          t(model_sol$ell1.T[1:model_sol$n.Z,drop=FALSE])%*%
                                          Z[[i-1]])),scale=param$mu_T)
    
    Z[[nb.simul+4]][i,]    <-
      rgamma(nb.traj,rpois(nb.traj,pmax(0,param$a_H/param$mu_H+
                                          t(model_sol$ell1.H[1:model_sol$n.Z,drop=FALSE])%*%
                                          Z[[i-1]])),scale=param$mu_H)
    
    eta[1:model_sol$n.eta,]<-apply(eta,2,function(i)param$Phi%*%i)+
      rnorm(model_sol$n.eta*nb.traj)
    W                      <-rbind(eta,
                                   Z[[nb.simul+1]][i,],
                                   Z[[nb.simul+2]][i,],
                                   Z[[nb.simul+3]][i,],
                                   Z[[nb.simul+4]][i,])
    
    #Climate variables
    Z_1.0 <-apply(Z[[i-1]],2,
                  function(x){A1[[i-1]]%*%x+omega0[[i-1]]})
    Wz    <-apply(W,2,function(x){omega[[i-1]]%*%x})
    Z[[i]]<-Z_1.0+Wz
    X[[i]]<-rbind(Z[[i]],W)
  }
  
  #extract path over time and remove time 0 (i.e., X_0)
  delc   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[1,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  y_tilde<-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[2,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  E      <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[3,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  E_ind  <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[4,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  Forc   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[5,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  M_at   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[6,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  M_up   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[7,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  M_lo   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[8,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  T_at   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[9,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  T_lo   <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[10,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  Cum_D  <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[11,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  Cum_E  <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[12,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  Cum_dc <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[13,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  H      <-matrix(t(extract(lapply(Z[1:nb.simul],function(x) x[14,]),1:nb.traj)),
                  ncol=nb.traj)[-1,]
  D      <- Z[[nb.simul+1]][-1,]
  N      <- Z[[nb.simul+2]][-1,]
  T_atW  <- Z[[nb.simul+3]][-1,]
  HW     <- Z[[nb.simul+4]][-1,]
  X      <-X[-1]
  
  vec_date<-seq(model_sol$vec_date[2],by=model_sol$tstep,length=nb.simul-1)
  
  
  mylist<-c(list("delc"=delc,"y_tilde"=y_tilde,"E"=E,"E_ind"=E_ind,"Forc"=Forc,
                 "M_at"=M_at,"M_up"=M_up,"M_lo"=M_lo,"T_at"=T_at,"T_lo"=T_lo,
                 "Cum_D"=Cum_D,"Cum_E"=Cum_E,"Cum_dc"=Cum_dc,"H"=H),
            list("D"=D,"N"=N,"T_atW"=T_atW,"HW"=HW),
            list("X"=X))
  return(mylist)
}


#* LAPLACE TRANSFORMS (LapT) FCTS ----
#CHANGE if add gamma0
#-------------------------------------Proposition: LapT of W
#*Functions for 'a.w.t' and 'b.w.t'
#*U is a dim(W)*N vector
#*

a1.w.fct <- function(model_sol,U,t){
  param <-model_sol$parameters
  n.eta <-model_sol$n.eta
  
  u.eta<-matrix(U[1:n.eta,],n.eta,dim(U)[2])
  u.d  <-U[n.eta+1,]
  u.n  <-U[n.eta+2,]
  u.T  <-U[n.eta+3,]
  u.H  <-U[n.eta+4,]
  
  a.w<-0
  a.w[1:dim(U)[2]]<- 0.5*matrix(apply(u.eta*u.eta,2,sum),
                                dim(U)[2],1)+
    matrix(u.d*param$a_D/(1-u.d*param$mu_D),
           dim(U)[2],1)+
    matrix((param$kappa_N^(t+1))*u.n*param$a_N/(1-u.n*param$mu_N),
           dim(U)[2],1)+
    matrix(0,
           dim(U)[2],1)+
    matrix(u.H*param$a_H/(1-u.H*param$mu_H),
           dim(U)[2],1)
  
  a.w<-matrix(a.w,dim(U)[2],1)
  
  if(!is.complex(u.d)|!is.complex(u.n)){
    a.w[u.d>=1/param$mu_D]     <- NaN
    a.w[u.n>=1/param$mu_N]     <- NaN
    a.w[u.T>=1/param$mu_T]     <- NaN
  }
  
  return(a.w)
}

b1.w.fct<-function(model_sol,U,t){
  param <-model_sol$parameters
  n.eta <-model_sol$n.eta
  
  u.eta<-matrix(U[1:n.eta,],n.eta,dim(U)[2])
  u.d  <-matrix(U[n.eta+1,],ncol=dim(U)[2])
  u.n  <-matrix(U[n.eta+2,],ncol=dim(U)[2])
  u.T  <-matrix(U[n.eta+3,],ncol=dim(U)[2])
  u.H  <-matrix(U[n.eta+4,],ncol=dim(U)[2])
  
  a                            <- matrix(0,
                                         model_sol$n.Z+model_sol$n.W,dim(U)[2])
  a[(model_sol$n.Z+1):
      (model_sol$n.Z+n.eta),]  <- t(param$Phi)%*%u.eta
  
  b                            <- matrix(model_sol$ell1.D%o%
                                           (u.d*param$mu_D/(1-u.d*param$mu_D))+
                                           model_sol$ell1.N%o%
                                           (u.n*param$mu_N*(param$kappa_N^(t+1))/
                                              (1-u.n*param$mu_N))+
                                           model_sol$ell1.T%o%
                                           (u.T*param$mu_T/(1-u.T*param$mu_T))+
                                           model_sol$ell1.H%o%
                                           (u.H*param$mu_H/(1-u.H*param$mu_H)),
                                         model_sol$n.Z+model_sol$n.W,
                                         dim(U)[2])
  
  b.w <- a + b
  
  if(!is.complex(u.d)|!is.complex(u.n)){
    b.w[,u.d>=1/param$mu_D]     <- NaN
    b.w[,u.n>=1/param$mu_N]     <- NaN
    b.w[,u.T>=1/param$mu_T]     <- NaN
  }
  
  return(b.w)
}


#-------------------------------------Proposition: One-Period ahead LapT of X
#*MIN=0, MAX=98, with 99=inf period
#*U: n.X*N matrix

a1.fct<-function(model_sol,U,t){
  omega.0   <-model_sol$omega0
  omega     <-model_sol$omega
  
  Uz<-matrix(U[1:model_sol$n.Z,],model_sol$n.Z,dim(U)[2])
  Uw<-matrix(U[(model_sol$n.Z+1):(model_sol$n.Z+model_sol$n.W),],
             model_sol$n.W,dim(U)[2])
  
  a <-t(Uz)%*%omega.0[[t+1]]+a1.w.fct(model_sol,Uw+t(omega[[t+1]])%*%Uz,t)
  
  
  return(a)
}

b1.fct<-function(model_sol,U,t){
  A.1       <-model_sol$A1
  omega     <-model_sol$omega
  
  Uz<-matrix(U[1:model_sol$n.Z,],model_sol$n.Z,dim(U)[2])
  Uw<-matrix(U[(model_sol$n.Z+1):model_sol$n.X,],model_sol$n.W,dim(U)[2])
  
  b <-rbind(t(A.1[[t+1]])%*%Uz,matrix(0,model_sol$n.W,dim(U)[2]))+
    b1.w.fct(model_sol,Uw+t(omega[[t+1]])%*%Uz,t)
  
  return(b)
}


#* INFINITE AND GN ----

#----------------------Functions for infinite a,b
#*U: matrix dim(U)=n.X*N
a1.w.fct.inf<-function(model_sol,U){
  param <-model_sol$parameters
  n.eta <-model_sol$n.eta
  
  u.eta<-matrix(U[1:n.eta,],n.eta,dim(U)[2])
  u.d  <-U[n.eta+1,]
  u.n  <-U[n.eta+2,]
  u.T  <-U[n.eta+3,]
  u.H  <-U[n.eta+4,]
  
  a.w<-0
  a.w[1:dim(U)[2]]<- 0.5*matrix(apply(u.eta*u.eta,2,sum),
                                dim(U)[2],1)+
    matrix(u.d*param$a_D/(1-u.d*param$mu_D),dim(U)[2],1)+
    matrix(0,dim(U)[2],1)+
    matrix(u.H*param$a_H/(1-u.H*param$mu_H),dim(U)[2],1)
  
  a.w<-matrix(a.w,dim(U)[2],1)
  
  if(!is.complex(u.d)|!is.complex(u.n)){
    a.w[u.d>=1/param$mu_D]     <- NaN
    a.w[u.n>=1/param$mu_N]     <- NaN
    a.w[u.T>=1/param$mu_T]     <- NaN
    a.w[u.H>=1/param$mu_H]     <- NaN
  }
  
  return(a.w)
}

b1.w.fct.inf<-function(model_sol,U){
  param <-model_sol$parameters
  n.eta <-model_sol$n.eta
  
  u.eta<-matrix(U[1:n.eta,],n.eta,dim(U)[2])
  u.d  <-matrix(U[n.eta+1,],ncol=dim(U)[2])
  u.n  <-matrix(U[n.eta+2,],ncol=dim(U)[2])
  u.T  <-matrix(U[n.eta+3,],ncol=dim(U)[2])
  u.H  <-matrix(U[n.eta+4,],ncol=dim(U)[2])
  
  a                            <- matrix(0,
                                         model_sol$n.Z+model_sol$n.W,dim(U)[2])
  a[(model_sol$n.Z+1):
      (model_sol$n.Z+n.eta),]  <- t(param$Phi)%*%u.eta
  
  b                            <- matrix(model_sol$ell1.D%o%
                                           (u.d*param$mu_D/(1-u.d*param$mu_D))+
                                           model_sol$ell1.T%o%
                                           (u.T*param$mu_T/(1-u.T*param$mu_T))+
                                           model_sol$ell1.H%o%
                                           (u.H*param$mu_H/(1-u.H*param$mu_H)),
                                         model_sol$n.Z+model_sol$n.W,
                                         dim(U)[2])
  
  b.w <- a + b
  
  if(!is.complex(u.d)|!is.complex(u.n)){
    b.w[,u.d>=1/param$mu_D]     <- NaN
    b.w[,u.n>=1/param$mu_N]     <- NaN
    b.w[,u.T>=1/param$mu_T]     <- NaN
  }
  
  return(b.w)
}

a1.fct.inf<-function(model_sol,U){
  omega.0.inf  <-model_sol$omega0.inf
  omega.inf    <-model_sol$omega.inf
  
  Uz<-matrix(U[1:model_sol$n.Z,],model_sol$n.Z,dim(U)[2])
  Uw<-matrix(U[(model_sol$n.Z+1):(model_sol$n.Z+model_sol$n.W),],model_sol$n.W,
             dim(U)[2])
  
  a.inf<-t(Uz)%*%omega.0.inf+a1.w.fct.inf(model_sol,Uw+t(omega.inf)%*%Uz)
  
  return(a.inf)
}

b1.fct.inf<-function(model_sol,U){
  omega.inf  <-model_sol$omega.inf
  A.1.inf    <-model_sol$A1.inf
  
  Uz<-matrix(U[1:model_sol$n.Z],model_sol$n.Z,1)
  Uw<-matrix(U[(model_sol$n.Z+1):(model_sol$n.Z+model_sol$n.W)],model_sol$n.W,1)
  
  b.inf<-rbind(t(A.1.inf)%*%Uz,matrix(0,model_sol$n.W,1))+
    b1.w.fct.inf(model_sol,Uw+t(omega.inf)%*%Uz)
  
  return(b.inf)
}


#-------------------------------------Gauss Newton function for mu_u1
#*x0: matrix of dimension "n.X*1", initial guess of infinite values
#*return nb of iterations, dev, mu_u1.inf, mu_c1, list dev

GN.function<-function(model_sol,x0){
  param     <-model_sol$parameters
  A1.inf    <-model_sol$A1.inf
  omega.inf <-model_sol$omega.inf
  #Jacobian
  e1.X    <-matrix(0,model_sol$n.Z+model_sol$n.W,1)
  e1.X[1] <-1
  mu_c1   <-e1.X
  mu_u1   <-x0
  eps     <-param$eps.GN
  J       <-matrix(0,model_sol$n.Z+model_sol$n.W,model_sol$n.Z+model_sol$n.W)
  dev     <-matrix(1,model_sol$n.Z+model_sol$n.W,1)
  ite     <-0
  listdev <-list(dev)
  tol     <-max(abs(dev*param$tol.GN))
  while((max(abs(dev))>tol)&(ite<50)){
    b     <-  b1.fct.inf(model_sol,(1-param$gamma)*(mu_u1+mu_c1))
    dev   <- -mu_u1+param$delta/(1-param$gamma)*b
    mu_u1 <-  param$delta/(1-param$gamma)*b
    ite           <-ite+1
    listdev[[ite]]<-dev
  }
  mylist<-list("ite"=ite,"dev"=dev,"mu_u1"=mu_u1,"listdev"=listdev,
               "mu_c1"=mu_c1)
  return(mylist)
}















### New pricing functions where varphi is evaluated in parallel
### at different omega's. For only one maturity at a time.


multi.lt.fct.Uh_MULTI<-function(model_sol,Uh,X=model_sol$X,t=0){
  H  <-dim(Uh)[3]
  K  <-dim(Uh)[2]
  U  <-list(matrix(Uh[,,H],ncol=K))
  if((t+H)>(model_sol$Tmax-1)){
    a.h<-matrix(rbind(c(a1.fct.inf(model_sol,U[[1]])),
                      matrix(NaN,H-1,K)),H,K)   
    i<-1
    #as long as t+H-i is larger than 98=max(a1/b1), we stay time-indep
    while((t+H-i)>(model_sol$Tmax-2)&i<H){
      U[[i+1]]<-Uh[,,H-i]+b1.fct.inf(model_sol,U[[i]])
      a.h[i+1]<-c(a1.fct.inf(model_sol,U[[i+1]]))
      i<-i+1
    }
    if(i<=(H-1)){
      for (k in i:(H-1)){
        U[[k+1]]<-Uh[,,H-k]+b1.fct(model_sol,U[[k]],t+H-k)
        a.h[k+1]<-c(a1.fct(model_sol,U[[k+1]],t+H-(k+1)))
      } 
    }
  }else{
    a.h<-matrix(rbind(c(a1.fct(model_sol,U[[1]],t+H-1)),
                      matrix(NaN,H-1,K)),H,K) 
    if(H>1){
      for (k in 1:(H-1)){
        U[[k+1]]<-Uh[,,H-k]+b1.fct(model_sol,U[[k]],t+H-k)
        a.h[k+1,]<- c(a1.fct(model_sol,U[[k+1]],t+H-(k+1)))
      }
    }
  }
  if(t>(model_sol$Tmax-2)){
    uX_t.h<-exp(apply(a.h,2,sum)+t(b1.fct.inf(model_sol,U[[H]]))%*%X)
    psi.0 <-apply(a.h,2,sum)
    psi.1 <-b1.fct.inf(model_sol,U[[H]])
  }else{
    uX_t.h<-exp(apply(a.h,2,sum)+t(b1.fct(model_sol,U[[H]],t))%*%X)
    psi.0 <-apply(a.h,2,sum)
    psi.1 <-b1.fct(model_sol,U[[H]],t) 
  }
  mylist<-list("psi.0"=psi.0,"psi.1"=psi.1,"uX_t.h"=uX_t.h)
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
  
  P.a.pi<-matrix(NaN,H-1,1)
  P.b.pi<-matrix(NaN,H-1,model_sol$n.X)
  
  U <- array(0,c(model_sol$n.X,K,H))
  
  if(H>1){
    for(k in 1:(H-1)){
      P.a.pi[k] <-a1.fct(model_sol,matrix(P.pi[t+k+1,],ncol=1),t+k)
      P.b.pi[k,]<-b1.fct(model_sol,matrix(P.pi[t+k+1,],ncol=1),t+k)
      
      U[,,k] <- - eta1[t+k+1,] - P.b.pi[k,] + P.pi[t+k,]
    }
  }
  U[,,H] <- P.pi[t+H,] + omega
  
  res.psi <- multi.lt.fct.Uh_MULTI(model_sol,Uh=U,X=X,t=t)
  
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
  #fx<-outer(x,b,function(r,c)Im(matrix(s1[,1],ncol=1)*exp(-1i*r*c))/r)*dx[,1]
  #varphi.hat<-varphi(model_sol,omega,H,X,t)[[3]][H]/2-1/pi*sum(fx)
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




