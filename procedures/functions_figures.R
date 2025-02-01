# ==============================================================================
# FIGURE 3. Damage function
# Figure_Damage_comparison.pdf
# make_figure_Damage_comparison.R
# + FIGURE III.1. Model calibration
#   Figure_Calibration.pdf
#   make_figure_calibration.R
#     i)
make_figure_calibration_Damages <- function(model_sol,
                                            damage.values = seq(0.995,.005,length.out=30),
                                            nb.values.interpol=1000,
                                            xD = exp(seq(-5,5,length.out = 2000)),   # to compute Riemann sum
                                            T.2100 = seq(1.5,5,length.out=20),
                                            vector.of.CI = c(0.5,0.8,0.9,0.95),
                                            main.title = "2100 Fractional damages",
                                            trace_lines = TRUE,
                                            lwd.CR = 2,
                                            indic_add_SLR = FALSE){
  
  P.col.line <- brocolors("crayons")["Teal Blue"]
  P.col <- adjustcolor( P.col.line, alpha.f = 0.15)
  
  damage.values.interpol <- seq(damage.values[1],
                                damage.values[length(damage.values)],
                                length.out=nb.values.interpol)
  gammaD <- - log(damage.values)
  gammaD.interpol <- - log(damage.values.interpol)
  
  #Matrix of temperatures in X
  all.cdf <- matrix(NaN,length(gammaD.interpol),length(T.2100))
  all.pdf <- matrix(NaN,length(gammaD.interpol)-1,length(T.2100))
  
  for(i in 1:length(T.2100)){
    aux <- Fourier.psi(model_sol,gammaD,xD,T0=model_sol$vector.ini$ini_Tat,
                       Tstar = T.2100[i],
                       tstar=model_sol$horiz.2100,
                       psi.LT=LT.CumD) # NB: by default, LT.CumD uses SLR damages
    
    tmp <- splinefun(x=gammaD, y=aux, method="hyman")
    
    fitted.cdf.values <- tmp(gammaD.interpol)
    fitted.pdf.values <- diff(fitted.cdf.values)
    
    all.cdf[,i] <- fitted.cdf.values
    all.pdf[,i] <- fitted.pdf.values
  }
  
  CI.P <- confidence_intervals_across_horizons(all.cdf,
                                               values.of.variable = 1 - damage.values.interpol,
                                               nb.values.variable = length(gammaD.interpol),
                                               vector.of.CI = vector.of.CI)
  all.CI.P  <- CI.P$all.CI
  
  indic.median.damage <- apply(all.cdf,2,function(x){which(x>.5)[1]})
  median.damage       <- 1 - damage.values.interpol[indic.median.damage]
  
  plot(T.2100,median.damage,type="l",ylim=c(0,.4),col="white",
       xlab="Temperature in 2100 (in \u00B0C)",
       ylab="Damages (fractions of 2100 output)",las=1,
       main=main.title)
  for(i in length(vector.of.CI):1){
    # P
    polygon(c(T.2100,rev(T.2100)),
            c(all.CI.P[1,,i],rev(all.CI.P[2,,i])),
            col=P.col,border = NA)
  }
  
  # Compute average damages:
  cond.mean.damage          <- NULL
  cond.mean.damage.with.SLR <- NULL
  for(i in 1:length(T.2100)){
    cond.mean.damage <- c(cond.mean.damage,
                          1 - LT.CumD(model_sol,u=-1,tstar=model_sol$horiz.2100,
                                      T0 = model_sol$vector.ini$ini_Tat,
                                      Tstar = T.2100[i],
                                      indic_add_SLR=FALSE))}
  lines(T.2100,cond.mean.damage,lwd=lwd.CR,lty=2,col="black")
  
  if(indic_add_SLR){
    for(i in 1:length(T.2100)){
      cond.mean.damage.with.SLR <- c(cond.mean.damage.with.SLR,
                                     1 - LT.CumD(model_sol,u=-1,tstar=model_sol$horiz.2100,
                                                 T0 = model_sol$vector.ini$ini_Tat,
                                                 Tstar = T.2100[i],
                                                 indic_add_SLR=TRUE))}
    lines(T.2100,cond.mean.damage.with.SLR,lwd=lwd.CR,lty=1,col="black")
  }
  
  if(trace_lines){
    points(c(2,4),
           c(1-model_sol$target_vector["ECumD2"],1-model_sol$target_vector["ECumD4"]),
           pch=15,col="red")
  }
  
  if(trace_lines){
    legend("topleft",
           legend=c("Mean, with SLR","Mean, w/o SLR","Targeted mean"),
           lty=c(1,2,NaN),
           pch=c(NaN,NaN,15),
           col=c("black","black","red"),
           #cex=1.5,
           lwd=2,bty = "n")
  }
  
  return(cond.mean.damage.with.SLR)
}

# ==============================================================================
# FIGURE 9. Effect of \mu on the conditional distribution x_t|y_t \sim \gamma0 (y_t/\mu, \mu)
# Figure_gamma0.pdf
# make_figure_gamma0_distri.R

#* FOURIER transform for permafrost-related emissions:
FOURIER<-function(model,x,gamma){
  # works only for i^th variable
  dx <- matrix(x-c(0,x[1:length(x)-1]),length(x),1)
  s1 <- matrix(PSI(model,1i*x),ncol=1)
  fx <- outer(x,gamma,function(r,c)Im(s1[,1]*exp(-1i*r*c))/r)*
    dx[,1]
  f <- 1/2-1/pi*apply(fx,2,sum)
  return(f)
}

PSI <- function(model,u){
  return(exp(model$y * u /(1 - u * model$mu) ))
}

# ==============================================================================
# FIGURE III.1. Model calibration
# Figure_Calibration.pdf
# make_figure_calibration.R
#  ii)
make_figure_calibration_N <- function(model_sol,
                                      CumN.values = seq(0,3000,length.out=60),
                                      nb.values.interpol=1000,
                                      xN = exp(seq(-10,10,length.out = 2000)),   # to compute Riemann sum
                                      T.2100 = seq(1.5,5,length.out=20),
                                      vector.of.CI = c(0.5,0.8,0.9,0.95),
                                      main.title = "2100 Permafrost-related emissions"){
  
  P.col.line <- brocolors("crayons")["Teal Blue"]
  P.col<-adjustcolor( P.col.line, alpha.f = 0.15)
  
  CumN.values.interpol <- seq(CumN.values[1],
                              CumN.values[length(CumN.values)],
                              length.out=nb.values.interpol)
  
  all.cdf <- matrix(NaN,length(CumN.values.interpol),length(T.2100))
  all.pdf <- matrix(NaN,length(CumN.values.interpol)-1,length(T.2100))
  
  for(i in 1:length(T.2100)){
    aux <- Fourier.psi(model_sol,CumN.values,xN,
                       T0=model_sol$vector.ini$ini_Tat,
                       Tstar = T.2100[i],
                       tstar=model_sol$horiz.2100,
                       psi.LT=LT.CumN)
    
    tmp <- splinefun(x=CumN.values, y=aux, method="hyman")
    
    fitted.cdf.values <- tmp(CumN.values.interpol)
    fitted.pdf.values <- diff(fitted.cdf.values)
    
    all.cdf[,i] <- fitted.cdf.values
    all.pdf[,i] <- fitted.pdf.values
  }
  
  CI.P <- confidence_intervals_across_horizons(all.cdf,
                                               values.of.variable = CumN.values.interpol,
                                               nb.values.variable = length(CumN.values.interpol),
                                               vector.of.CI = vector.of.CI)
  all.CI.P  <- CI.P$all.CI
  
  indic.median.CumN <- apply(all.cdf,2,function(x){which(x>.5)[1]})
  median.CumN <- CumN.values.interpol[indic.median.CumN]
  
  plot(T.2100,median.CumN,type="l",
       ylim=c(0,model_sol$target_vector["ECumN4"]+3*model_sol$target_vector["stdCumN4"]),
       col="white",
       xlab="Temperature in 2100 (in \u00B0C)",
       ylab="Permafrost-related emissions (in GtC)",las=1,
       main=main.title)
  
  for(i in length(vector.of.CI):1){
    # P
    polygon(c(T.2100,rev(T.2100)),
            c(all.CI.P[1,,i],rev(all.CI.P[2,,i])),
            col=P.col,border = NA)
  }
  
  res <- Param.Gamma0.CumN(model_sol,T0=model_sol$vector.ini$ini_Tat,
                           Tstar=T.2100,tstar=model_sol$horiz.2100)
  cond.mean.CumN <- res$lambda * res$mu
  lines(T.2100,cond.mean.CumN,lwd=2,lty=2)
  
  points(c(2,4),
         c(model_sol$target_vector["ECumN2"],model_sol$target_vector["ECumN4"]),
         pch=15,col="red")
  
  legend("topleft",
         legend=c("Mean","Targeted mean"),
         lty=c(2,NaN),
         pch=c(NaN,NaN,15),
         col=c("black","red"),
         #cex=1.5,
         lwd=2,bty = "n")
  
  return(1) 
}

# iii)
make_figure_calibration_SLR <- function(model_sol,
                                        H.values = seq(0,3,length.out=50),
                                        nb.H = 1000,
                                        xH = exp(seq(-5,5,length.out = 2000)),  # to compute Riemann sum
                                        T.2100 = seq(1.5,5,length.out=20),
                                        vector.of.CI = c(0.5,0.8,0.9,0.95),
                                        main.title = "2100 Sea-level rise"){
  
  t.star <- model_sol$horiz.2100
  
  P.col.line <- brocolors("crayons")["Teal Blue"]
  P.col<-adjustcolor( P.col.line, alpha.f = 0.15)
  
  H.values.interpol <- seq(H.values[1],
                           H.values[length(H.values)],
                           length.out=nb.H)
  
  all.cdf <- matrix(NaN,nb.H,length(T.2100))
  all.pdf <- matrix(NaN,nb.H-1,length(T.2100))
  
  for(i in 1:length(T.2100)){
    aux <- Fourier.psi(model_sol,
                       H.values,xH,
                       T0=model_sol$vector.ini$ini_Tat,
                       Tstar = T.2100[i],
                       tstar=model_sol$horiz.2100,
                       psi.LT=LT.SLR)
    
    tmp <- splinefun(x=H.values, y=aux, method="hyman")
    
    fitted.cdf.values <- tmp(H.values.interpol)
    fitted.pdf.values <- diff(fitted.cdf.values)
    
    all.cdf[,i] <- fitted.cdf.values
    all.pdf[,i] <- fitted.pdf.values
  }
  
  CI.P <- confidence_intervals_across_horizons(all.cdf,
                                               values.of.variable = H.values.interpol,
                                               nb.values.variable = length(H.values.interpol),
                                               vector.of.CI = vector.of.CI)
  all.CI.P  <- CI.P$all.CI
  
  indic.median.H <- apply(all.cdf,2,function(x){which(x>.5)[1]})
  median.H       <- H.values.interpol[indic.median.H]
  
  plot(T.2100,median.H,type="l",
       ylim=range(H.values),col="white",
       xlab="Temperature in 2100 (in \u00B0C)",
       ylab="SLR (in meters)",las=1,
       main=main.title)
  
  for(i in length(vector.of.CI):1){
    # P
    polygon(c(T.2100,rev(T.2100)),
            c(all.CI.P[1,,i],rev(all.CI.P[2,,i])),
            col=P.col,border = NA)
  }
  
  # Add median:
  indic.median.H <- apply(all.cdf,2,function(x){which(x>.5)[1]})
  median.H <- H.values.interpol[indic.median.H]

  # Add expectation:
  res <- Param.Gamma0.H(model_sol,
                        T0=model_sol$vector.ini$ini_Tat,
                        Tstar=T.2100,
                        tstar=model_sol$horiz.2100)
  cond.mean.H <- res$lambda * res$mu
  lines(T.2100,cond.mean.H,lwd=2,lty=2)
  
  points(c(2,4),c(model_sol$target_vector[["EH2"]],
                  model_sol$target_vector[["EH4"]]),pch=15,col="red")
  
  legend("topleft",
         legend=c("Mean","Targeted mean"),
         lty=c(2,NaN),
         pch=c(NaN,15),
         col=c("black","red"),
         #cex=1.5,
         lwd=2,bty = "n")
  
  return(1)
}

# ==============================================================================
# FIGURE 6. The term structure of real rates
# Figure_YC_RF.pdf
# + FIGURE V.3. Comparison of temperature risk premiums
#   Figure_TRP_comparison.pdf
# make_figure_YC_RF.R
# ==============================================================================

#   i)
compute_CMT <- function(model_sol,EV,Price.ZC,
                        maturities = c(2,10), # in number of periods
                        nb = model_sol$horiz.2100+1,
                        theta0 = list(c(log(0.17),-1/21*log(0.17)))
){
  
  h_cst.1 <- maturities[1]
  h_cst.2 <- maturities[2]
  
  # Model with climate change (CC)
  cc.1 <-compute_cst_h(model_sol,h_cst.1,nb)
  cc.2 <-compute_cst_h(model_sol,h_cst.2,nb)
  
  # Model without climate change risk but agent does not re-optimize
  model_ncc <- model_sol
  model_ncc$parameters$mu_D <- 0.00001
  model_ncc$parameters$mu_H <- 0.00001
  model_ncc$parameters$mu_N <- 0.00001
  model_ncc$parameters$mu_T <- 0.00001
  # and agents re-optimize:
  model_ncc <- model_solve(model_ncc)
  ncc.1 <-compute_cst_h(model_ncc,h_cst.1,nb)
  ncc.2 <-compute_cst_h(model_ncc,h_cst.2,nb)
  
  # Model without climate change risk + mu_c=E_t(delc)
  model_ncc.mu_c     <- model_ncc
  model_ncc.mu_c$parameters$a_D  <- 0
  model_ncc.mu_c$parameters$b_D  <- 0
  model_ncc.mu_c$parameters$b_sk <- 0
  E_delc <- EV.fct(model_sol,h=model_sol$Tmax-1)$EX$delc
  model_ncc.mu_c$mu_c <- c(model_sol$X[1],E_delc)
  model_ncc.mu_c      <- model_solve(model_ncc.mu_c,indic_mitig = FALSE)
  
  ncc.mu_c.1 <-compute_cst_h(model_ncc.mu_c,h_cst.1,nb)
  ncc.mu_c.2 <-compute_cst_h(model_ncc.mu_c,h_cst.2,nb)
  
  return(list(cc.1=cc.1,cc.2=cc.2,
              ncc.1=ncc.1,ncc.2=ncc.2,
              ncc.mu_c.1=ncc.mu_c.1,ncc.mu_c.2=ncc.mu_c.2))
} 

#  ii)
make_figure_CMT <- function(expected.yds,
                            model_sol,
                            maturities = c(2,10),
                            indic_only_first=FALSE){
  
  cc.1       <- expected.yds$cc.1
  cc.2       <- expected.yds$cc.2
  ncc.1      <- expected.yds$ncc.1
  ncc.2      <- expected.yds$ncc.2
  ncc.mu_c.1 <- expected.yds$ncc.mu_c.1
  ncc.mu_c.2 <- expected.yds$ncc.mu_c.2
  
  # colors for 3 models
  col<-c("black",
         "honeydew4",
         "dodgerblue3")
  
  ylim=c(min(cc.1[,2],ncc.1[,2],ncc.mu_c.1[,2],cc.2[,2],ncc.2[,2],ncc.mu_c.2[,2])-.25,
         max(cc.1[,2],ncc.1[,2],ncc.mu_c.1[,2],cc.2[,2],ncc.2[,2],ncc.mu_c.2[,2])+.25)
  
  #Plot
  cexaxs  <- 1
  cexlab  <- 1
  cexmain <- 1.2
  
  plot(cc.1[,1],cc.1[,2],
       type="b", xlab="Date",ylab="Interest rate (in percent)",
       ylim=ylim,
       las=1,
       cex.main=cexmain,cex.axis=cexaxs,cex.lab=cexlab,
       main=paste(ifelse(indic_only_first,"(b)","(a)")," Expected path of ",
                  model_sol$tstep*maturities[1],"-year real interest rate",sep=""),
       lwd=2,col=col[1],pch=1)
  lines(ncc.1[,1],ncc.1[,2],type="b",
        lwd=2,col=col[2],pch=2)
  lines(ncc.mu_c.1[,1],ncc.mu_c.1[,2],type="b",
        lwd=2,col=col[3],pch=4)
  grid()
  
  legend("bottomleft",
         legend=c("Baseline",
                  "Deterministic damages",
                  "Deterministic damages, same growth path as baseline"),
         lty=rep(1,4),
         pch=c(1,2,4),
         pt.bg = 'white',
         col=col,
         lwd=rep(2,4),
         bty = "n",cex=1)
  
  if(!indic_only_first){
    plot(cc.2[,1],cc.2[,2],
         type="b", xlab="Date",ylab="Interest rate (in percent)",
         ylim=ylim,
         las=1,
         cex.main=cexmain,cex.axis=cexaxs,cex.lab=cexlab,
         main=paste("(b) Expected path of ",model_sol$tstep*maturities[2],"-year real interest rate",sep=""),
         lwd=2,col=col[1],pch=1)
    lines(ncc.2[,1],ncc.2[,2],type="b",
          lwd=2,col=col[2],pch=2)
    lines(ncc.mu_c.2[,1],ncc.mu_c.2[,2],type="b",
          lwd=2,col=col[3],pch=4)
    grid()
  }
  
  return(1) 
}

