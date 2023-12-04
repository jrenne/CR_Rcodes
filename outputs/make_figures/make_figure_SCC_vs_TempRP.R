# ==============================================================================
# Relationship between climate beta and temperature risk premium
# ==============================================================================


all.SCC.RP <- NULL
all.T.RP   <- NULL

omega_ZCB  <- matrix(0,model_sol$n.X)
omega_T.at <- matrix(0,model_sol$n.X)
omega_T.at[which(model_sol$names.var.X=="T_at")] <- 1
omega_H <- matrix(0,model_sol$n.X)
omega_H[which(model_sol$names.var.X=="H")] <- 1

all.multip.factor <- seq(0.01,1.5,length.out=7)
indic.baseline    <- which((all.multip.factor - 1)^2==min((all.multip.factor - 1)^2))

cl <- makeCluster(number.of.cores)
registerDoParallel(cl)

save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))

clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))


FILE = paste("/outputs/Figures/Figure_SCCvsTRP.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=7, height=3)

par(mfrow=c(1,3))
par(plt=c(.2,.95,.2,.85))

cases <- seq(0,1,by=.5)
i.case <- 0
for(indic.CRRA in cases){
  
  i.case <- i.case + 1
  print(paste("   * Case ",i.case," (out of ",length(cases),")",sep=""))
  
  if(indic.CRRA==1){
    gamma <- 1.45
    xlim  <- c(-.025,.02)
    ylim  <- c(-.07,.03)
  }else if(indic.CRRA == .5){
    gamma <- 1.001
    xlim  <- c(-.025,.02)
    ylim  <- c(-.07,.03)
  }else{
    gamma <- model_sol$parameters$gamma
    xlim  <- c(-.1,.3)
    ylim  <- c(-1,.4)
  }
  
  if((indic.CRRA==1)){
    eval(parse(text = gsub(" "," ",
                           paste("main.t <- expression(paste('(c) Power utility, ',gamma,' = ',",
                                 gamma,",sep=''))",sep="")
    )))
  }else{
    if((indic.CRRA==.5)){
      eval(parse(text = gsub(" "," ",
                             paste("main.t <- expression(paste('(b) Epstein-Zin (~CRRA), ',gamma,' = ',",
                                   gamma,",sep=''))",sep="")
      )))
    }else{
      eval(parse(text = gsub(" "," ",
                             paste("main.t <- expression(paste('(a) Epstein-Zin, ',gamma,' = ',",
                                   gamma,",sep=''))",sep="")
      )))
    }
  }
  
  if(indic.CRRA==0){
    ylab <- "SCC minus NPV of benefits (as a fraction of SCC)"
  }else{
    ylab <- ""
  }
  
  plot(0,0,col="white",
       xlab="2100 Temperature risk premium, in °C",
       ylab=ylab,
       main=main.t,
       xlim=xlim,
       ylim=ylim)

  for(indic.variable in c("conso vol","damages")){
    
    all.res <- foreach(i = 1:length(all.multip.factor),.combine=rbind) %dopar% {
      
      model_new <- model_sol
      model_new$parameters$gamma <- gamma
      
      targets <- model_sol$target_vector
      
      multip.factor <- all.multip.factor[i]
      
      if(indic.variable=="conso vol"){
        targets["sigma_c0"] <- multip.factor * model_sol$target_vector["sigma_c0"]
      }
      if(indic.variable=="damages"){
        targets["ECumD2"] <- 1 - multip.factor * 
          (1 - model_sol$target_vector["ECumD2"])
        targets["ECumD4"] <- 1 - multip.factor * 
          (1 - model_sol$target_vector["ECumD4"])
        targets["stdCumD4"] <- min(multip.factor,1) * model_sol$target_vector["stdCumD4"]
        model_new$parameters$b_sk <- multip.factor * model_sol$parameters$b_sk
      }
      
      model_new$target_vector <- targets
      model_new <- solveParam4D(model_new)
      model_new <- solveParam4H(model_new)
      model_new <- solveParam4N(model_new)
      model_new <- solveParam4c(model_new,
                                indic_CRRA = (indic.CRRA==1))
      
      model_sol_new <- model_solve(model_new,
                                   indic_CRRA = (indic.CRRA==1))
      
      if(!(indic.CRRA==1)){
        SCC <- scc.fct(model_sol_new,h=0)
      }else{
        SCC <- scc.fct.CRRA(model_sol_new)$SCC.CO2
      }
      
      EV    <- EV.fct(model_sol_new,h=H)
      
      ZCB <- varphi(model_sol_new,omega_ZCB,H)
      
      ET.P  <- EV$EX$T_at[1:H]
      ET.Q  <- varphi.tilde(model_sol_new,omega_T.at,H)[[1]]/ZCB$P.t
      
      EH.P  <- EV$EX$H[1:H]
      EH.Q  <- varphi.tilde(model_sol_new,omega_H,H)[[1]]/ZCB$P.t
      
      T.RP <- c(ET.Q - ET.P)
      H.RP <- c(EH.Q - EH.P)
      
      # computation of climate beta
      
      HH <- 200 # maximum horizon considered in infinite sums
      epsilon <- .1
      X.shock <- model_sol_new$X
      X.shock[which(model_sol$names.var.X=="M_at")] <- - epsilon + 
        X.shock[which(model_sol$names.var.X=="M_at")]
      # prices.C   <- varphi(model_sol_new,
      #                      omega.varphi = omega_C,
      #                      H = 100)
      # prices.C.shock   <- varphi(model_sol_new,
      #                            omega.varphi = omega_C,
      #                            H = 100,
      #                            X = X.shock)
      # D <- sum(prices.C$P.t) - sum(prices.C.shock$P.t)
      
      EC       <- NULL
      EC.shock <- NULL
      for(h in 1:HH){
        Uh <- matrix(model_sol_new$mu_c1,model_sol$n.X,h)
        res.lt       <- multi.lt.fct.Uh(model_sol_new,Uh,X=model_sol_new$X,t=0)
        res.lt.shock <- multi.lt.fct.Uh(model_sol_new,Uh,X=X.shock,t=0)
        EC <- c(EC,res.lt$uX_t.h)
        EC.shock <- c(EC.shock,res.lt.shock$uX_t.h)
      }
      # Discounted expected benefits:
      ZCB <- varphi(model_sol_new,omega_ZCB,HH)
      NPV.CO2 <- sum(ZCB$P.t * (EC.shock - EC)/epsilon) *
        10^3 * model_sol$parameters$c0 / 3.667
      
      c(SCC,NPV.CO2,ET.Q[H],ET.P[H])
    }
    
    SCC.RP <- (all.res[,1] - all.res[,2])/all.res[,1]
    ET.RP  <- all.res[,3] - all.res[,4]
    
    
    if(indic.variable=="conso vol"){
      col <- "black"
      lty <- 1
      
      all.SCC.RP <- rbind(all.SCC.RP,SCC.RP)
      all.T.RP   <- rbind(all.T.RP,ET.RP)
    }
    if(indic.variable=="damages"){
      col <- "dark grey"
      lty <- 3
    }
    
    abline(h=0,col="grey",lty=3)
    abline(v=0,col="grey",lty=3)
    lines(ET.RP,SCC.RP,lwd=2,col=col,lty=lty)
    
    points(ET.RP[1],SCC.RP[1],pch=1,lwd=2,cex=1.5,col=col)
    points(ET.RP[indic.baseline],
           SCC.RP[indic.baseline],pch=0,lwd=2,cex=1.5,col=col)
    points(ET.RP[length(ET.RP)],
           SCC.RP[length(ET.RP)],pch=3,lwd=2,cex=1.5,col=col)
    
    if(indic.CRRA==0){
      if(indic.variable=="conso vol"){
        legend("bottomright",
               legend=c(expression(paste(sigma[a]," = 0",sep="")),
                        expression(paste("Baseline value of ",sigma[A],sep="")),
                        expression(paste("Largest value of ",sigma[A]," (x 1.5)",sep=""))),
               lty=c(NaN,NaN),
               col=col,
               pch=c(1,0,3),
               cex=1.2,
               lwd=c(2,2,2),bty = "n")
      }
      if(indic.variable=="damages"){
        legend("topright",
               legend=c("No damages","Baseline damages","Large damages (x 1.5)"),
               lty=c(NaN,NaN),
               col=col,
               pch=c(1,0,3),
               cex=1.2,
               lwd=c(2,2,2),bty = "n")
      }
    }
  }
}

dev.off()

# To be saved:
all.SCC.RP.CR <- all.SCC.RP
all.T.RP.CR   <- all.T.RP
save(all.SCC.RP.CR,all.T.RP.CR,
     file="outputs/results/SCC_vs_TRP_CR.Rdat")


stopCluster(cl)
file.remove("outputs/toto.Rdata")

