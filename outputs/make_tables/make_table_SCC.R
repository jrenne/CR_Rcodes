# ==============================================================================
# TABLE 5. Social Cost of Carbon, in US $ per ton of CO2
# table_SCC_sensitiv.txt
# TABLE S.1. 2100 Temperature risk premium
# table_TRP_sensitiv.txt
# TABLE S.2. Long-term interest rate (maturity: 2100)
# table_LTR_sensitiv.txt
# ==============================================================================

# Define format of figures:
nb.dec <- 2 # number of decimal numbers
format.nb  <- paste("%.",nb.dec,"f",sep="")
format.nb0 <- paste("%.",0,"f",sep="")
format.nb1 <- paste("%.",1,"f",sep="")
format.nb2 <- paste("%.",2,"f",sep="")
format.nb3 <- paste("%.",3,"f",sep="")
format.nb4 <- paste("%.",4,"f",sep="")
format.nb5 <- paste("%.",5,"f",sep="")
make.entry <- function(x,format.nb,dollar=TRUE){
  if(dollar){
    output <- paste("$",sprintf(format.nb,x),"$",sep="")
  }else{
    output <- paste(sprintf(format.nb,x),sep="")
  }
  return(output)
}

omega_ZCB <- matrix(0,model_sol$n.X,1)
omega_T.at <- omega_ZCB
omega_H    <- omega_ZCB
omega_T.at[which(model_sol$names.var.X=="T_at")] <- 1
omega_H[which(model_sol$names.var.X=="H")]       <- 1

H <- model_sol$horiz.2100
matur <- H


targets <- model_sol$target_vector
targets["delta"] <- NaN
targets["nu"]    <- NaN

# Define new models:

targets_smallD <- targets
targets_smallD["ECumD2"]   <- 1 - (1-targets["ECumD2"])/2
targets_smallD["ECumD4"]   <- 1 - (1-targets["ECumD4"])/2
targets_smallD["stdCumD4"] <- targets["stdCumD4"]/2
targets_smallD["EH2"]      <- targets["EH2"]/2
targets_smallD["EH4"]      <- targets["EH4"]/2
targets_smallD["stdH4"]    <- targets["stdH4"]/2

targets_noD.uncert <- targets
targets_noD.uncert["stdCumD4"] <- .00001

targets_smallN <- targets
targets_smallN["ECumN2"]   <- .00001
targets_smallN["ECumN4"]   <- .00001
targets_smallN["ECumNinf"] <- .0001
targets_smallN["stdCumN4"] <- sqrt(.00001)

targets_noN.uncert <- targets
targets_noN.uncert["stdCumN4"] <- .00001

targets_noT.uncert <- targets
targets_noT.uncert["stdTat2100"] <- NaN

targets_smallH <- targets
targets_smallH["EH2"] <- .00001
targets_smallH["EH4"] <- .00002

targets_noH.uncert <- targets
targets_noH.uncert["stdH4"] <- .00001

targets_smallGrowth0 <- targets
targets_smallGrowth0["mu_c0"] <- targets["mu_c0"]*4

targets_noGrowth0.uncert <- targets
targets_noGrowth0.uncert["sigma_c0"] <- .00001

targets_no.uncert <- targets
targets_no.uncert["stdCumD4"]   <- .00001
targets_no.uncert["stdCumN4"]   <- .00001
targets_no.uncert["stdH4"]      <- .00001
targets_no.uncert["sigma_c0"]   <- .00001
targets_no.uncert["stdTat2100"] <- NaN

targets_smallD_no.uncert <- targets
targets_smallD_no.uncert["stdCumD4"]   <- .000001
targets_smallD_no.uncert["stdCumN4"]   <- .000001
targets_smallD_no.uncert["stdH4"]      <- .000001
targets_smallD_no.uncert["sigma_c0"]   <- .000001
targets_smallD_no.uncert["stdTat2100"] <- NaN
targets_smallD_no.uncert["ECumD2"]     <- 1 - (1-targets["ECumD2"])/2
targets_smallD_no.uncert["ECumD4"]     <- 1 - (1-targets["ECumD4"])/2
targets_smallD_no.uncert["EH2"]        <- targets["EH2"]/2
targets_smallD_no.uncert["EH4"]        <- targets["EH4"]/2
targets_smallD_no.uncert["ECumN2"]     <- .00001
targets_smallD_no.uncert["ECumN4"]     <- .00001
targets_smallD_no.uncert["ECumNinf"]   <- .0001
targets_smallD_no.uncert["stdCumN4"]   <- sqrt(.00001)

targets_delta_15 <- targets
targets_delta_15["delta"] <- (1- .015)^5

targets_delta_05 <- targets
targets_delta_05["delta"] <- (1- .005)^5

targets_nu_2 <- targets
targets_nu_2["nu"] <- 2

targets_nu_4 <- targets
targets_nu_4["nu"] <- 4

all.targets <- rbind(targets,
                     targets_delta_15,
                     targets_delta_05,
                     targets_nu_4,
                     targets_nu_2,
                     targets_no.uncert,
                     targets_noGrowth0.uncert,
                     targets_noD.uncert,
                     targets_noN.uncert,
                     targets_noH.uncert,
                     targets_noT.uncert,
                     targets_smallD,
                     targets_smallN,
                     targets_smallH
)
rownames(all.targets) <- c("Baseline",
                           "Large rate of pref. for present ($-\\log(\\delta)/\\Delta t=1.5\\%$)",
                           "Small rate of pref. for present ($-\\log(\\delta)/\\Delta t=0.5\\%$)",
                           "Large climate sensitivity ($\\nu=4$)",
                           "Small climate sensitivity ($\\nu=2$)",
                           "Deterministic model",
                           "No exog. techno. uncertainty ($\\sigma_A=0$)",
                           "No uncert. on damage specif. ($\\mu_D=0$)",
                           "No uncert. on permafrost specif. ($\\mu_N=0$)",
                           "No uncert. on sea-level specif. ($\\mu_H=0$)",
                           "No uncert. on temperature specif. ($\\mu_T=0$)",
                           "Small damages (twice lower)",
                           "No permaf. releases ($a^{(N)}=b^{(N)}=0$)",
                           "No sea-level-specific damages ($a^{(H)}=b^{(H)}=0$)"
)


for(indic.CRRA in c(FALSE,TRUE)){# if FALSE, use Epstein-Zin, otherwise CRRA
  if(indic.CRRA){
    values.of.gamma <- c(1.5,1.001,2)
    print("   * Time-separable CRRA preferences (about 2 minutes)")
  }else{
    print("   * Epstein-Zin preferences (about 1 minute)")
    values.of.gamma <- c(model_sol$parameters$gamma,2,10)
  }
  
  all.specif <- cbind(values.of.gamma%x%matrix(1,dim(all.targets)[1],1),
                      matrix(1,length(values.of.gamma),1) %x% all.targets)
  
  # Run parallel computations:
  cl <- makeCluster(number.of.cores)
  registerDoParallel(cl)
  save.image("outputs/toto.Rdata")
  clusterEvalQ(cl,load("outputs/toto.Rdata"))
  clusterEvalQ(cl,library(MASS))
  clusterEvalQ(cl,library(expm))
  all_results <- foreach(i = 1:dim(all.specif)[1], 
                         .combine=rbind) %dopar% {
                           
                           targets <- all.specif[i,2:dim(all.specif)[2]]
                           names(targets) <- c(names(model_sol$target_vector),
                                               "delta","nu")
                           
                           gamma <- all.specif[i,1]
                           
                           model_new <- model_sol
                           
                           model_new$parameters$gamma <- gamma
                           
                           model_new$target_vector <- targets
                           
                           model_new <- solveParam4D(model_new)
                           model_new <- solveParam4H(model_new)
                           model_new <- solveParam4N(model_new)
                           model_new <- solveParam4c(model_new,
                                                     indic_CRRA = indic.CRRA)
                           if(is.na(model_new$target_vector["stdTat2100"])){
                             model_new$parameters$mu_T <- .000001
                           }
                           if(!is.na(model_new$target_vector["delta"])){
                             model_new$parameters$delta <- model_new$target_vector["delta"]
                           }
                           if(!is.na(model_new$target_vector["nu"])){
                             model_new$parameters$nu <- model_new$target_vector["nu"]
                           }
                           
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
                           X.shock[which(model_sol$names.var.X=="M_at")] <- - epsilon + 
                             X.shock[which(model_sol$names.var.X=="M_at")]
                           
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
                           
                           c(SCC,NPV.CO2,T.RP,H.RP,ZCB$r.t)
                         }
  
  stopCluster(cl)
  file.remove("outputs/toto.Rdata")
  
  all_SCC <- all_results[,1]
  all_SCC <- matrix(all_SCC,dim(all.targets)[1],length(values.of.gamma))
  rownames(all_SCC) <- rownames(all.targets)
  colnames(all_SCC) <- paste("gamma = ",values.of.gamma,sep="")
  #print(all_SCC)
  
  all_NPV <- all_results[,2]
  all_BETA <- 100*(all_SCC - all_NPV)/all_SCC
  all_BETA <- matrix(all_BETA,dim(all.targets)[1],length(values.of.gamma))
  rownames(all_BETA) <- rownames(all.targets)
  colnames(all_BETA) <- paste("gamma = ",values.of.gamma,sep="")
  #print(all_BETA)
  
  all_TRP <- all_results[,1+1+H]
  all_TRP <- matrix(all_TRP,dim(all.targets)[1],length(values.of.gamma))
  rownames(all_TRP) <- rownames(all.targets)
  colnames(all_TRP) <- paste("gamma = ",values.of.gamma,sep="")
  #print(all_TRP)
  
  all_HRP <- all_results[,1+1+2*H]
  all_HRP <- matrix(all_HRP,dim(all.targets)[1],length(values.of.gamma))
  rownames(all_HRP) <- rownames(all.targets)
  colnames(all_HRP) <- paste("gamma = ",values.of.gamma,sep="")
  #print(all_HRP)
  
  all_LTR <- all_results[,1+1+2*H+matur]
  all_LTR <- matrix(all_LTR,dim(all.targets)[1],length(values.of.gamma))
  rownames(all_LTR) <- rownames(all.targets)
  colnames(all_LTR) <- paste("gamma = ",values.of.gamma,sep="")
  #print(all_LTR)
  
  # Save results: --------------------------------------------------------------
  save(all_SCC,all_LTR,all_HRP,all_TRP,all_NPV,all_BETA,
       all.targets,values.of.gamma,
       file=paste("outputs/results/sensitivity_SCC_",
                  ifelse(indic.CRRA,"CRRA","EZ"),
                  ".Rdat",sep=""))
}





# Prepare Latex tables: --------------------------------------------------------


# SCC --------------------------------------------------------------------------

latex.table <- NULL

for(indic.CRRA in c(FALSE,TRUE)){
  
  # load results: --------------------------------------------------------------
  load(file=paste("outputs/results/sensitivity_SCC_",
                  ifelse(indic.CRRA,"CRRA","EZ"),
                  ".Rdat",sep=""))
  
  if(!indic.CRRA){
    latex.table <- rbind(latex.table,
                         "\\hline",
                         "\\multicolumn{7}{l}{{\\qquad\\bf A. Epstein-Zin preferences}}\\\\",
                         "\\hline")
  }else{
    latex.table <- rbind(latex.table,
                         "\\hline",
                         "\\multicolumn{7}{l}{{\\qquad\\bf B. Time-separable CRRA preferences}}\\\\",
                         "\\hline")
  }
  
  this.line <- "Specification"
  for(j in 1:length(values.of.gamma)){
    this.line <- paste(this.line,"&\\multicolumn{2}{c}{$\\gamma=",
                       values.of.gamma[j],"$}",sep="")
  }
  latex.table <- rbind(latex.table,
                       paste(this.line,"\\\\",sep=""),
                       "\\hline")
  
  for(i in 1:dim(all.targets)[1]){
    this.line <- rownames(all.targets)[i]
    for(j in 1:length(values.of.gamma)){
      if(i==1){
        this.line <- paste(this.line,"&{\\bf ",
                           round(all_SCC[i,j]),"}&",sep="")
      }else{
        rel.diff <- (all_SCC[i,j] - all_SCC[1,j])/all_SCC[1,j]
        this.line <- paste(this.line,"&",
                           round(all_SCC[i,j]),
                           "&$",ifelse(rel.diff*sign(all_SCC[1,j])>=0,"+","-"),"$\\textit{",make.entry(100*abs(rel.diff),format.nb1,dollar = FALSE),"\\%}",
                           sep="")
      }
    }
    latex.table <- rbind(
      latex.table,
      paste(this.line,"\\\\",sep=""))
  }
}

# Save results:
latex.file <- paste("outputs/Tables/table_SCC_sensitiv.txt",sep="")
write(latex.table, file = latex.file)




# TRP --------------------------------------------------------------------------

latex.table <- NULL

for(indic.CRRA in c(FALSE,TRUE)){
  
  # load results: --------------------------------------------------------------
  load(file=paste("outputs/results/sensitivity_SCC_",
                  ifelse(indic.CRRA,"CRRA","EZ"),
                  ".Rdat",sep=""))
  
  if(!indic.CRRA){
    latex.table <- rbind(latex.table,
                         "\\hline",
                         "\\multicolumn{7}{l}{{\\qquad\\bf A. Epstein-Zin preferences}}\\\\",
                         "\\hline")
  }else{
    latex.table <- rbind(latex.table,
                         "\\hline",
                         "\\multicolumn{7}{l}{{\\qquad\\bf B. Time-separable CRRA preferences}}\\\\",
                         "\\hline")
  }
  
  this.line <- "Specification"
  for(j in 1:length(values.of.gamma)){
    this.line <- paste(this.line,"&\\multicolumn{2}{c}{$\\gamma=",
                       values.of.gamma[j],"$}",sep="")
  }
  latex.table <- rbind(latex.table,
                       paste(this.line,"\\\\",sep=""),
                       "\\hline")
  
  for(i in 1:dim(all.targets)[1]){
    this.line <- rownames(all.targets)[i]
    for(j in 1:length(values.of.gamma)){
      if(i==1){
        this.line <- paste(this.line,"&{\\bf ",
                           make.entry(all_TRP[i,j],format.nb3,dollar = FALSE),"}&",sep="")
      }else{
        rel.diff <- all_TRP[i,j] - all_TRP[1,j]
        this.line <- paste(this.line,"&",
                           make.entry(all_TRP[i,j],format.nb3),
                           "&$",ifelse(rel.diff>=0,"+","-"),"$\\textit{",make.entry(abs(rel.diff),format.nb3,dollar = FALSE),"\\degree C}",
                           sep="")
      }
    }
    latex.table <- rbind(
      latex.table,
      paste(this.line,"\\\\",sep=""))
  }
}
# Save results:
latex.file <- paste("outputs/Tables/table_TRP_sensitiv.txt",sep="")
write(latex.table, file = latex.file)




# LTR --------------------------------------------------------------------------

latex.table <- NULL

for(indic.CRRA in c(FALSE,TRUE)){
  
  # load results: --------------------------------------------------------------
  load(file=paste("outputs/results/sensitivity_SCC_",
                  ifelse(indic.CRRA,"CRRA","EZ"),
                  ".Rdat",sep=""))
  
  if(!indic.CRRA){
    latex.table <- rbind(latex.table,
                         "\\hline",
                         "\\multicolumn{7}{l}{{\\qquad\\bf A. Epstein-Zin preferences}}\\\\",
                         "\\hline")
  }else{
    latex.table <- rbind(latex.table,
                         "\\hline",
                         "\\multicolumn{7}{l}{{\\qquad\\bf B. Time-separable CRRA preferences}}\\\\",
                         "\\hline")
  }
  
  this.line <- "Specification"
  for(j in 1:length(values.of.gamma)){
    this.line <- paste(this.line,"&\\multicolumn{2}{c}{$\\gamma=",
                       values.of.gamma[j],"$}",sep="")
  }
  latex.table <- rbind(latex.table,
                       paste(this.line,"\\\\",sep=""),
                       "\\hline")
  
  for(i in 1:dim(all.targets)[1]){
    this.line <- rownames(all.targets)[i]
    for(j in 1:length(values.of.gamma)){
      if(i==1){
        this.line <- paste(this.line,"&{\\bf ",
                           make.entry(all_LTR[i,j],format.nb2,dollar = FALSE),"}&",sep="")
      }else{
        rel.diff <- all_LTR[i,j] - all_LTR[1,j]
        this.line <- paste(this.line,"&",
                           make.entry(all_LTR[i,j],format.nb2,dollar = FALSE),
                           "&$",ifelse(rel.diff>=0,"+","-"),"$\\textit{",round(100*abs(rel.diff))," bps}",
                           sep="")
      }
    }
    latex.table <- rbind(
      latex.table,
      paste(this.line,"\\\\",sep=""))
  }
}
# Save results:
latex.file <- paste("outputs/Tables/table_LTR_sensitiv.txt",sep="")
write(latex.table, file = latex.file)