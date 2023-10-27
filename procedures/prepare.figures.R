
make_figure_calibration_Damages <- function(model_sol,
                                            damage.values = seq(0.995,.005,length.out=30),
                                            nb.values.interpol=1000,
                                            xD = exp(seq(-5,5,length.out = 2000)),   # to compute Riemann sum
                                            T.2100 = seq(1.5,5,length.out=20),
                                            vector.of.CI = c(0.5,0.8,0.9,0.95),
                                            main.title = "2100 Fractional damages"){
  
  P.col.line <- brocolors("crayons")["Teal Blue"]
  P.col<-adjustcolor( P.col.line, alpha.f = 0.15)
  
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
                       psi.LT=LT.CumD)
    
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
       ylab="Fractional damage",las=1,
       main=main.title)
  for(i in length(vector.of.CI):1){
    # P
    polygon(c(T.2100,rev(T.2100)),
            c(all.CI.P[1,,i],rev(all.CI.P[2,,i])),
            col=P.col,border = NA)
  }
  
  lines(T.2100,median.damage,lwd=2,lty=2)
  
  # Compute average damages:
  cond.mean.damage <- NULL
  for(i in 1:length(T.2100)){
    cond.mean.damage <- c(cond.mean.damage,
                          1 - LT.CumD(model_sol,u=-1,tstar=model_sol$horiz.2100,
                                      T0=model_sol$vector.ini$ini_Tat,Tstar = T.2100[i]))}
  lines(T.2100,cond.mean.damage,lwd=2)
  
  points(c(2,4),
         c(1-model_sol$target_vector["ECumD2"],1-model_sol$target_vector["ECumD4"]),
         pch=15,col="red")
  
  legend("topleft",
         legend=c("mean","median","target (mean)"),
         lty=c(1,2,NaN),
         pch=c(NaN,NaN,15),
         col=c("black","black","red"),
         #cex=1.5,
         lwd=2,bty = "n")
  
  return(1)
}




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
  
  plot(T.2100,median.CumN,type="l",ylim=c(0,model_sol$target_vector["ECumN4"]+3*model_sol$target_vector["stdCumN4"]),col="white",
       xlab="Temperature in 2100 (in \u00B0C)",
       ylab="Permafrost-related emissions (in GtC)",las=1,
       main=main.title)
  
  for(i in length(vector.of.CI):1){
    # P
    polygon(c(T.2100,rev(T.2100)),
            c(all.CI.P[1,,i],rev(all.CI.P[2,,i])),
            col=P.col,border = NA)
  }
  
  
  lines(T.2100,median.CumN,lwd=2,lty=2)
  
  res <- Param.Gamma0.CumN(model_sol,T0=model_sol$vector.ini$ini_Tat,
                           Tstar=T.2100,tstar=model_sol$horiz.2100)
  cond.mean.CumN <- res$lambda * res$mu
  lines(T.2100,cond.mean.CumN,lwd=2)
  
  points(c(2,4),
         c(model_sol$target_vector["ECumN2"],model_sol$target_vector["ECumN4"]),
         pch=15,col="red")
  
  legend("topleft",
         legend=c("mean","median","target (mean)"),
         lty=c(1,2,NaN),
         pch=c(NaN,NaN,15),
         col=c("black","black","red"),
         #cex=1.5,
         lwd=2,bty = "n")
  
  return(1) 
}



make_figure_calibration_SLR <- function(model_sol,
                                        H.values = seq(0,3,length.out=50),
                                        nb.H = 1000,
                                        xH = exp(seq(-5,5,length.out = 2000)),   # to compute Riemann sum
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
  lines(T.2100,median.H,lwd=2,lty=2)
  
  # Add expectation:
  H.2100 <- model_sol$vector.ini$ini_H + 
    model_sol$parameters$a_H * t.star +
    model_sol$parameters$b_H * (.5*((t.star-1)*T.2100+
                                      (t.star+1)*model_sol$vector.ini$ini_Tat))
  lines(T.2100,H.2100,lwd=2)
  
  points(c(2,4),c(model_sol$target_vector[["EH2"]],
                  model_sol$target_vector[["EH4"]]),pch=15,col="red")
  
  legend("topleft",
         legend=c("mean","median","target (mean)"),
         lty=c(1,2,NaN),
         pch=c(NaN,NaN,15),
         col=c("black","black","red"),
         #cex=1.5,
         lwd=2,bty = "n")
  
  return(1)
}





make_figure_trajectory_and_pdf <- function(model_sol,EV,Price.ZC,
                                           variable,# example: "T_at"
                                           year.density, #example: 2100
                                           x = exp(seq(-5,5,length.out = 1000)),
                                           values, # where Fourier is performed
                                           name.of.variable, # title for chart
                                           unit # for chart's label
){
  
  H <- length(EV$EX[[1]])
  H.specif <- which(EV$date==year.density)
  indic.variable <- which(variable==model_sol$names.var.X)
  
  # Prepare omega vectors (for pricing):
  omega_ZCB <- matrix(0,model_sol$n.X)
  omega_Variable <- omega_ZCB
  omega_Variable[which(model_sol$names.var.X==variable)] <- 1
  
  # Compute expected trajectories (under P and Q):
  ET.P    <- EV$EX[[indic.variable]][1:H]
  ET.Q    <- varphi.tilde(model_sol,omega_Variable,H)[[1]]/Price.ZC
  
  #distribution and expected temperature
  
  cdf.Q<-c(varphi.hat.fast(model_sol,omega = omega_ZCB,
                           H=H.specif,x,a = omega_Variable,
                           b=values)/Price.ZC[H.specif])
  cdf.P <- fourier(model_sol,x,values,H.specif,indic.variable)
  
  cdf.P <- pmax(pmin(cdf.P,1),0)
  cdf.Q <- pmax(pmin(cdf.Q,1),0)
  for(i in 2:length(cdf.P)){
    if(cdf.P[i]<cdf.P[i-1]){
      cdf.P[i] <- cdf.P[i-1]
    }
    if(cdf.Q[i]<cdf.Q[i-1]){
      cdf.Q[i] <- cdf.Q[i-1]
    }
  }
  
  length.extended <- 500
  values.extended <- seq(values[1],values[length(values)],
                         length.out=length.extended)
  tmp.P               <- splinefun(x=values, y=cdf.P, method="hyman")
  fitted.cdf.P.values <- tmp.P(values.extended)
  fitted.pdf.P.values <- diff(fitted.cdf.P.values)
  tmp.Q               <- splinefun(x=values, y=cdf.Q, method="hyman")
  fitted.cdf.Q.values <- tmp.Q(values.extended)
  fitted.pdf.Q.values <- diff(fitted.cdf.Q.values)
  
  P.col.line <- "#18a7b5"
  Q.col.line <- "#ff8243"
  P.col<-adjustcolor( P.col.line, alpha.f = 0.15)
  Q.col<-adjustcolor( Q.col.line, alpha.f = 0.15)
  
  ylim=c(min(ET.P,ET.Q) - .2*(max(ET.P,ET.Q) - min(ET.P,ET.Q)),
         max(ET.P,ET.Q) + .6*(max(ET.P,ET.Q) - min(ET.P,ET.Q)))
  ylim <- c(values[1],values[length(values)])
  
  par(plt=c(.15,1,.25,.85))
  plot(EV$date,ET.P,col=P.col.line,lwd=2,type="l",ylim=ylim,
       xlab="horizon",ylab=unit,main=name.of.variable,
       cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
  lines(EV$date,ET.Q,col=Q.col.line,lwd=2)
  abline(v=year.density,lty=3)
  abline(h=ET.P[H.specif],lty=3,col=P.col.line)
  abline(h=ET.Q[H.specif],lty=3,col=Q.col.line)
  
  par(plt=c(0,.95,.25,.85))
  plot(fitted.pdf.P.values,values.extended[2:length.extended],
       type="l",lwd=2,col=P.col.line,ylim=ylim,xaxt="n",yaxt="n",xlab="",ylab="",
       cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
  lines(fitted.pdf.Q.values,values.extended[2:length.extended],
        lwd=2,col=Q.col.line)
  # abline(h=values.extended[which.min(abs(fitted.cdf.Q.values-0.5))],
  #        col=Q.col.line,lty=1,lwd=1)                                              #median of Q
  abline(h=ET.Q[H.specif],col=Q.col.line,lty=3,lwd=1)                                    #mean of Q
  # abline(h=values.extended[which.min(abs(fitted.cdf.P.values-0.5))],
  #        col=P.col.line,lty=1,lwd=1)                                              #median of P
  abline(h=ET.P[H.specif],col=P.col.line,lty=3,lwd=1)                                    #mean of P
  abline(v=0,col="grey")
  
  polygon(c(fitted.pdf.P.values,rev(fitted.pdf.P.values)),
          c(values.extended[2:length.extended],0*values.extended[2:length.extended]),
          col=P.col,border = NULL)
  polygon(c(fitted.pdf.Q.values,rev(fitted.pdf.Q.values)),
          c(values.extended[2:length.extended],0*values.extended[2:length.extended]),
          col=Q.col,border = NULL)
  
  legend("topright",
         legend=c("Physical pdf for selected date","Risk-adjusted pdf for selected date"),
         lty=c(1,1),
         col=c(P.col.line,Q.col.line),cex=1.5,
         lwd=c(2,2),bty = "n")
  
  return(1)
}



# H <- 30
# EV <- EV.fct(model_sol,h=H)
# # Compute P and Q proba:
# Price.ZC <- varphi(model_sol,omega_ZCB,H = H)[[3]]
# 
# 
# par(mfrow=c(2,2))
# 
# make_figure_trajectory_and_pdf(model_sol,EV,Price.ZC,
#                                variable = "T_at",
#                                year.density = 2100,
#                                values = seq(1,EV$EX$T_at[H]+5*sqrt(EV$VX[[TAT]][H]),by=.1), # where Fourier is performed
#                                name.of.variable = expression(paste("Atmospheric temperature",T[AT],sep="")), # title for chart
#                                unit = "in °C")
# 
# make_figure_trajectory_and_pdf(model_sol,EV,Price.ZC,
#                                variable = "H",
#                                year.density = 2100, 
#                                values = seq(0,2.5,by=.05), # where Fourier is performed
#                                name.of.variable = expression(paste("Global mean sea level ",H,sep="")), # title for chart
#                                unit = "in meters")
# 




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
  model_ncc$parameters$mu_D <- 0.0001
  model_ncc$parameters$mu_H <- 0.0001
  model_ncc$parameters$mu_N <- 0.0001
  model_ncc$parameters$mu_T <- 0.0001
  # and agents re-optimize:
  model_ncc <- model_solve(model_ncc,theta0)
  ncc.1 <-compute_cst_h(model_ncc,h_cst.1,nb)
  ncc.2 <-compute_cst_h(model_ncc,h_cst.2,nb)
  
  # Model without climate change risk + mu_c=E_t(delc)
  model_ncc.mu_c     <- model_ncc
  model_ncc.mu_c$parameters$a_D  <- 0
  model_ncc.mu_c$parameters$b_D  <- 0
  model_ncc.mu_c$parameters$b_sk <- 0
  E_delc <- EV.fct(model_sol,h=model_sol$Tmax-1)$EX$delc
  model_ncc.mu_c$mu_c<- c(model_sol$X[1],E_delc)
  model_ncc.mu_c     <- model_solve(model_ncc.mu_c,theta0,indic_mitig = FALSE)
  
  ncc.mu_c.1 <-compute_cst_h(model_ncc.mu_c,h_cst.1,nb)
  ncc.mu_c.2 <-compute_cst_h(model_ncc.mu_c,h_cst.2,nb)
  
  return(list(cc.1=cc.1,cc.2=cc.2,
              ncc.1=ncc.1,ncc.2=ncc.2,
              ncc.mu_c.1=ncc.mu_c.1,ncc.mu_c.2=ncc.mu_c.2))
} 


make_figure_CMT <- function(expected.yds,
                            model_sol,
                            maturities = c(2,10)){
  
  cc.1       <- expected.yds$cc.1
  cc.2       <- expected.yds$cc.2
  ncc.1      <- expected.yds$ncc.1
  ncc.2      <- expected.yds$ncc.2
  ncc.mu_c.1 <- expected.yds$ncc.mu_c.1
  ncc.mu_c.2 <- expected.yds$ncc.mu_c.2
  
  # colors for 3 models
  col<-c(brocolors("crayons")["Orange Red"],
         brocolors("crayons")["Wild Blue Yonder"],
         brocolors("crayons")["Violet Blue"])
  
  ylim=c(min(cc.1[,2],ncc.1[,2],ncc.mu_c.1[,2],cc.2[,2],ncc.2[,2],ncc.mu_c.2[,2])-.5,
         max(cc.1[,2],ncc.1[,2],ncc.mu_c.1[,2],cc.2[,2],ncc.2[,2],ncc.mu_c.2[,2])+.5)
  
  #Plot
  cexaxs  <- 1
  cexlab  <- 1
  cexmain <- 1
  
  plot(cc.1[,1],cc.1[,2],
       type="b", xlab="Date",ylab="Interest rate (in percent)",
       ylim=ylim,
       las=1,
       cex.main=cexmain,cex.axis=cexaxs,cex.lab=cexlab,
       main=paste("(a) ",model_sol$tstep*maturities[1],"-year real interest rate",sep=""),
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
  
  plot(cc.2[,1],cc.2[,2],
       type="b", xlab="Date",ylab="Interest rate (in percent)",
       ylim=ylim,
       las=1,
       cex.main=cexmain,cex.axis=cexaxs,cex.lab=cexlab,
       main=paste("(b) ",model_sol$tstep*maturities[2],"-year real interest rate",sep=""),
       lwd=2,col=col[1],pch=1)
  lines(ncc.2[,1],ncc.2[,2],type="b",
        lwd=2,col=col[2],pch=2)
  lines(ncc.mu_c.2[,1],ncc.mu_c.2[,2],type="b",
        lwd=2,col=col[3],pch=4)
  grid()
  
  return(1) 
}

