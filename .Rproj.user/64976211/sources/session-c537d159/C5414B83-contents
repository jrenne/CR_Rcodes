
# Do something with mu_T (at least in no uncertainty at all)

values.of.gamma <- c(model$parameters$gamma,
                     model$parameters$gamma - 3,
                     model$parameters$gamma + 3)

targets <- model_sol$target_vector

targets_smallD <- targets
targets_smallD["ECumD2"] <- 1 - (1-targets["ECumD2"])/2
targets_smallD["ECumD4"] <- 1 - (1-targets["ECumD4"])/2
targets_smallD["stdCumD4"] <- targets["stdCumD4"]/2

targets_noD.uncert <- targets
targets_noD.uncert["stdCumD4"] <- .00001

targets_smallN <- targets
targets_smallN["ECumN2"] <- .00001
targets_smallN["ECumN4"] <- .00001
targets_smallN["stdCumN4"] <- sqrt(.00001)

targets_noN.uncert <- targets
targets_noN.uncert["stdCumN4"] <- .00001

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
targets_smallD_no.uncert["stdCumD4"]   <- .00001
targets_smallD_no.uncert["stdCumN4"]   <- .00001
targets_smallD_no.uncert["stdH4"]      <- .00001
targets_smallD_no.uncert["sigma_c0"]   <- .00001
targets_smallD_no.uncert["stdTat2100"] <- NaN
targets_smallD_no.uncert["ECumD2"] <- 1 - (1-targets["ECumD2"])/2
targets_smallD_no.uncert["ECumD4"] <- 1 - (1-targets["ECumD4"])/2

all.targets <- rbind(
  targets,
  targets_no.uncert,
  targets_noD.uncert,
  targets_smallD,
  targets_smallD_no.uncert,
  targets_noN.uncert,
  targets_smallN,
  #targets_noGrowth0.uncert,
  #targets_smallGrowth0,
  targets_noH.uncert,
  targets_smallH
)
rownames(all.targets) <- c("Baseline",
                           "Deterministic model",
                           "No damage uncertainty ($\\mu_T=0$)",
                           "Small damages (twice lower in expectation)",
                           "Small damages and deterministic",
                           "No permafrost uncertainty ($\\mu_N=0$)",
                           "No permaf.-based carbon releases ($a_N=b_N=0$)",
                           #"no growth uncertainty",
                           #"small growth",
                           "No sea-level uncertainty ($\\mu_H=0$)",
                           "No sea-level risks ($a_H=b_H=0$)")

all.specif <- cbind(values.of.gamma%x%matrix(1,dim(all.targets)[1],1),
                    matrix(1,length(values.of.gamma),1) %x% all.targets)

# Run parallel computations:
cl <- makeCluster(number.of.cores)
registerDoParallel(cl)
save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))
clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))
all_SCC <- foreach(i = 1:dim(all.specif)[1], 
                   .combine=rbind) %dopar% {
                     
                     targets <- all.specif[i,2:dim(all.specif)[2]]
                     names(targets) <- names(model$target_vector)
                     
                     gamma   <- all.specif[i,1]
                     
                     model_new <- model
                     
                     model_new$parameters$gamma <- gamma
                     
                     model_new$target_vector <- targets
                     model_new <- solveParam4D(model_new)
                     model_new <- solveParam4H(model_new)
                     model_new <- solveParam4N(model_new)
                     model_new <- solveParam4c(model_new)
                     if(is.na(model_new$target_vector["stdTat2100"])){
                       model_new$parameters$mu_T <- .00001
                     }
                     
                     model_sol_new <- model_solve(model_new,theta0)
                     
                     SCC <- scc.fct(model_sol_new,h=0)
                     
                     SCC
                   }

stopCluster(cl)

all_SCC <- matrix(all_SCC,dim(all.targets)[1],length(values.of.gamma))
rownames(all_SCC) <- rownames(all.targets)
colnames(all_SCC) <- paste("gamma = ",values.of.gamma,sep="")

print(all_SCC)


# Prepare Latex table:

this.line <- "Specification"
for(j in 1:length(values.of.gamma)){
  this.line <- paste(this.line,"&\\multicolumn{2}{c}{$\\gamma=",
                     values.of.gamma[j],"$}",sep="")
}
latex.table <- paste(this.line,"\\\\",sep="")
latex.table <- rbind(latex.table,
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
                         "&$",ifelse(rel.diff>=0,"+","-"),"$\\textit{",round(100*abs(rel.diff)),"\\%}",
                         sep="")
    }
  }
  latex.table <- rbind(
    latex.table,
    paste(this.line,"\\\\",sep=""))
}

latex.file <- paste("outputs/Tables/table_SCC_sensitiv.txt",sep="")
write(latex.table, file = latex.file)



