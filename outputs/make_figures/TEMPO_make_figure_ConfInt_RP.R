# ==============================================================================
# Determine main determinants of SCC
# Approach: moment targets and other key parameters are modified. In the end, we
#           run a regression where the SCC is regressed on targets & parameters
# ==============================================================================

max.chge <- .5 # maximum absolute change in target

Nb.draws <- 300

H <- model_sol$horiz.2100

# gammas <- c(model$parameters$gamma,2)
# deltas <- c(model$parameters$delta,(1-.01)^model$tstep)
# 
# rhos <- 1 - deltas^(1/model$tstep)
# 
# all.gamma <- gammas %x% c(1,1)
# all.delta <- c(1,1) %x% deltas

# Prepare omega vectors (for pricing):
omega_ZCB  <- matrix(0,model_sol$n.X)
omega_T.at <- matrix(0,model_sol$n.X)
omega_T.at[which(model_sol$names.var.X=="T_at")] <- 1
omega_H <- matrix(0,model_sol$n.X)
omega_H[which(model_sol$names.var.X=="H")] <- 1


cl <- makeCluster(number.of.cores)
registerDoParallel(cl)

save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))

clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))

all.res <- foreach(i = 1:Nb.draws, .combine=rbind) %dopar% {
  
  set.seed(i)
  
  model_sol_new <- model_sol
  
  # Delta:
  rho <- 1 - model_sol$parameters$delta^(1/model$tstep)
  rho_new <- rho * (1 + max.chge*runif(1,min=-1,max=1))
  model_sol_new$parameters$delta <- (1 - rho_new)^model$tstep
  
  # Gamma:
  model_sol_new$parameters$gamma <-
    model_sol_new$parameters$gamma * (1 + max.chge*runif(1,min=-1,max=1))
  
  # Damage targets:
  model_sol_new$target_vector["ECumD2"] <- 1 - (1-model_sol$target_vector["ECumD2"]) *
    (1 + max.chge*runif(1,min=-1,max=1))
  model_sol_new$target_vector["ECumD4"] <- model_sol_new$target_vector["ECumD2"] -
    (model_sol$target_vector["ECumD2"] - model_sol$target_vector["ECumD4"]) *
    (1 + max.chge*runif(1,min=-1,max=1))
  model_sol_new$target_vector["stdCumD4"] <- model_sol$target_vector["stdCumD4"] *
    (1 + max.chge*runif(1,min=-1,max=1))
  
  # Permafrost-related targets:
  model_sol_new$target_vector["ECumN2"] <- model_sol$target_vector["ECumN2"] *
    (1 + max.chge*runif(1,min=-1,max=1))
  model_sol_new$target_vector["ECumN4"] <- model_sol_new$target_vector["ECumN2"] +
    (model_sol$target_vector["ECumN4"] - model_sol$target_vector["ECumN2"]) *
    (1 + max.chge*runif(1,min=-1,max=1))
  model_sol_new$target_vector["stdCumN4"] <- model_sol$target_vector["stdCumN4"] *
    (1 + max.chge*runif(1,min=-1,max=1))
  model_sol_new$target_vector["ECumNinf"] <- model_sol$target_vector["ECumNinf"] *
    (1 + max.chge*runif(1,min=-1,max=1))
  
  # SLR targets:
  model_sol_new$target_vector["EH2"] <- model_sol$target_vector["EH2"] *
    (1 + max.chge*runif(1,min=-1,max=1))
  model_sol_new$target_vector["EH4"] <- model_sol_new$target_vector["EH2"] +
    (model_sol$target_vector["EH4"] - model_sol$target_vector["EH2"]) *
    (1 + max.chge*runif(1,min=-1,max=1))
  model_sol_new$target_vector["stdH4"] <- model_sol$target_vector["stdH4"] *
    (1 + max.chge*runif(1,min=-1,max=1))
  
  # Consumption moment targets:
  model_sol_new$target_vector["mu_c0"] <- model_sol$target_vector["mu_c0"] *
    (1 + max.chge*runif(1,min=-1,max=1))
  model_sol_new$target_vector["sigma_c0"] <- model_sol$target_vector["sigma_c0"] *
    (1 + max.chge*runif(1,min=-1,max=1))
  
  # Temperature dynamics uncertainty:
  mu_T <- model_sol_new$parameters$mu_T * (1 + max.chge*runif(1,min=-1,max=1))
  model_sol_new$parameters$mu_T <- mu_T
  
  targets_and_mu <- c(model_sol_new$target_vector,mu_T)
  save(targets_and_mu,
       file=paste("outputs/results/targets_and_mu",i,".Rdat",sep=""))
  
  model_sol_new <- solveParam4D(model_sol_new)
  model_sol_new <- solveParam4H(model_sol_new)
  model_sol_new <- solveParam4N(model_sol_new)
  model_sol_new <- solveParam4c(model_sol_new)
  
  model_sol_new <- model_solve(model_sol_new,indic_CRRA = FALSE)
  
  SCC <- scc.fct(model_sol_new,h=0)
  
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
  
  c(SCC,T.RP,H.RP,ZCB$r.t,
    model_sol_new$target_vector,mu_T,model_sol_new$parameters$gamma,rho_new)
}


colnames(all.res) <- c("scc",
                       paste("T.RP.",1:H,sep=""),
                       paste("H.RP.",1:H,sep=""),
                       paste("yld.",1:H,sep=""),
                       names(model$target_vector),
                       "mu_T","gamma","rho")

stopCluster(cl)
file.remove("outputs/toto.Rdata")

all.res <- all.res[complete.cases(all.res[,1]),]

scc <- all.res[,1]
X   <- all.res[,(1+3*H+1):dim(all.res)[2]]

eq <- lm(scc~X)
summary(eq)

X_red <- X[,c("ECumD2","ECumD4","stdCumD4","ECumN4",
              "gamma","rho")]
eq_red <- lm(scc~X_red)
summary(eq_red)





