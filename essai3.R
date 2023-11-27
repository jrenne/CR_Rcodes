
indic.CRRA <- 0

gamma <- 7

all.sigma.c0 <- seq(0,2*model_sol$target_vector["sigma_c0"],length.out=8)


tic()
cl <- makeCluster(number.of.cores)
registerDoParallel(cl)

save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))

clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))

all.res <- foreach(i = 1:length(all.sigma.c0), .combine=rbind) %dopar% {
  
  targets <- model_sol$target_vector
  
  #targets["stdCumD4"] <- model_sol$target_vector["stdCumD4"]*1
  #targets["ECumD2"] <- .975
  #targets["ECumD4"] <- .95
  
  targets["sigma_c0"] <- all.sigma.c0[i]
  
  model_new <- model
  model_new$parameters$gamma <- gamma
  
  #model_new$parameters$mu_T <- model$parameters$mu_T*2
  
  model_new$target_vector <- targets
  model_new <- solveParam4D(model_new)
  model_new <- solveParam4H(model_new)
  model_new <- solveParam4N(model_new)
  model_new <- solveParam4c(model_new,
                            indic_CRRA = indic.CRRA)
  
  model_sol_new <- model_solve(model_new,
                               indic_CRRA = indic.CRRA)
  
  if(!indic.CRRA){
    SCC <- scc.fct(model_sol_new,h=0)
  }else{
    SCC <- scc.fct.CRRA(model_sol_new)$SCC.CO2
  }
  
  EV    <- EV.fct(model_sol_new,h=H)
  
  ZCB <- varphi(model_sol_new,omega_ZCB,H)
  
  ET.P  <- EV$EX$T_at[1:H]
  ET.Q  <- 
    varphi.tilde(model_sol_new,omega_T.at,H)[[1]]/ZCB$P.t
  
  EH.P  <- EV$EX$H[1:H]
  EH.Q  <- 
    varphi.tilde(model_sol_new,omega_H,H)[[1]]/ZCB$P.t
  
  T.RP <- c(ET.Q - ET.P)
  H.RP <- c(EH.Q - EH.P)
  
  # computation of climate beta
  
  HH <- 200 # maximum horizon considered in infinite sums
  epsilon <- .1
  X.shock <- model_sol_new$X
  X.shock[which(model$names.var.X=="M_at")] <- - epsilon + 
    X.shock[which(model$names.var.X=="M_at")]
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
    10^3 * model$parameters$c0 / 3.667
  
  
  c(SCC,NPV.CO2,ET.Q[H],ET.P[H])
}

stopCluster(cl)
file.remove("outputs/toto.Rdata")
toc()

SCC.RP <- all.res[,1] - all.res[,2]
ET.RP  <- all.res[,3] - all.res[,4]

plot(ET.RP,SCC.RP,type="l",lwd=2)
abline(h=0,col="grey",lty=3)
abline(v=0,col="grey",lty=3)

