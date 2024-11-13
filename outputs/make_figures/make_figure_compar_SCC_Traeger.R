# ==============================================================================
# Figure comparing SCC with ACE model (Traeger, 2023)
# ==============================================================================

Damage_Traeger_4 <- 1 - compute_alternative_damage(T=4, type="L. Traeger")
Damage_CR_4 <- 1 - model$target_vector["ECumD4"]
factor_mult_Damage <- Damage_Traeger_4/Damage_CR_4

omega_ZCB <- matrix(0,model_sol$n.X)

all_rho <- as.numeric(levels(as.factor(data_ACE$prpp)))
all_scc_EZ   <- NULL
all_scc_CRRA <- NULL
all_scc_EZ_lowDamage   <- NULL
all_scc_CRRA_lowDamage <- NULL
all_5yrRates_EZ   <- NULL
all_5yrRates_CRRA <- NULL
all_5yrRates_EZ_lowDamage   <- NULL
all_5yrRates_CRRA_lowDamage <- NULL
for(lowDamage in c(TRUE,FALSE)){
  
  if(lowDamage){
    print("--- Low damages ---")
  }else{
    print("--- Baseline damages ---")
  }
  
  model_base <- model
  
  if(lowDamage){
    model_base$parameters$a_D  <- model_base$parameters$a_D * factor_mult_Damage
    model_base$parameters$b_D  <- model_base$parameters$b_D * factor_mult_Damage
    model_base$parameters$b_sk <- 0
  }
  
  for(rho in all_rho){
    print(rho)
    
    model_new <- model_base
    #model_new$target_vector["mu_c0"] <- .2
    model_new$parameters$delta <- (1 - rho)^model$tstep
    #model_new <- solveParam4c(model_new)
    #model_new <- solveParam4T(model_new)
    model_new_sol <- model_solve(model_new,indic_CRRA = FALSE)
    
    gamma <- 1.45
    model.CRRA <- model_base
    model.CRRA$parameters$gamma <- gamma
    model.CRRA$parameters$delta <- (1 - rho)^model$tstep
    model.CRRA <- solveParam4c(model.CRRA,indic_CRRA=TRUE)
    model_CRRA_sol <- model_solve(model.CRRA,
                                  indic_mitig = TRUE,
                                  indic_CRRA = TRUE)
    
    if(!lowDamage){
      all_scc_EZ <- c(all_scc_EZ,
                      scc.fct(model_new_sol,0))
      all_scc_CRRA <- c(all_scc_CRRA,
                        scc.fct.CRRA(model_CRRA_sol)$SCC.CO2)
      
      Price.ZC <- varphi(model_new_sol,
                         omega.varphi=omega_ZCB,
                         H = 1)
      all_5yrRates_EZ <- c(all_5yrRates_EZ,Price.ZC$r.t/100)
      Price.ZC <- varphi(model_CRRA_sol,
                         omega.varphi=omega_ZCB,
                         H = 1)
      all_5yrRates_CRRA <- c(all_5yrRates_CRRA,Price.ZC$r.t/100)
      
    }else{
      all_scc_EZ_lowDamage <- c(all_scc_EZ_lowDamage,
                                scc.fct(model_new_sol,0))
      all_scc_CRRA_lowDamage <- c(all_scc_CRRA_lowDamage,
                                  scc.fct.CRRA(model_CRRA_sol)$SCC.CO2)
      
      Price.ZC <- varphi(model_new_sol,
                         omega.varphi=omega_ZCB,
                         H = 1)
      all_5yrRates_EZ_lowDamage <- c(all_5yrRates_EZ_lowDamage,Price.ZC$r.t/100)
      Price.ZC <- varphi(model_CRRA_sol,
                         omega.varphi=omega_ZCB,
                         H = 1)
      all_5yrRates_CRRA_lowDamage <- c(all_5yrRates_CRRA_lowDamage,Price.ZC$r.t/100)
    }
  }
}


# Plot ----
FILE = "/outputs/Figures/Figure_SCC_ACE_comparison.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=11, width=6, height=4)

par(plt=c(.19,.95,.2,.95))

data_ACE <- read.csv("data/ACE_Traeger_2023.csv")

plot(log(all_scc_EZ),pch=15,col="red",ylim=c(2,12),
     cex=1.6,
     xlab=expression(paste("Pure rate of preference for present ",(1-delta)^{1/5},sep="")),
     ylab="",
     xaxt = "n", yaxt = "n",
     xlim=c(0.3,3.7))

title(ylab=expression(paste("SCC, in dollars per ton of ",CO[2],sep="")), line=5,
      cex.lab=1)


points(log(all_scc_CRRA),pch=17,col="red",cex=1.6)

points(log(all_scc_EZ_lowDamage),pch=0,col="red",cex=1.6)
points(log(all_scc_CRRA_lowDamage),pch=2,col="red",cex=1.6)


type_prpp <- apply(matrix(data_ACE$prpp,ncol=1),1,function(x){which(x==all_rho)})
points(type_prpp,log(data_ACE$SCC),pch=19,col="#00000033",cex=1.3)

# X-axis
axis(1, at = 1:3,labels = all_rho)
# Y-axis
ylabels <- log(c(0,1,10,25,100,250,1000,2500,10000,25000))
axis(2, at = ylabels, labels = exp(ylabels),las=1)
for(i in 1:length(ylabels)){
  abline(h=ylabels[i],col="grey",lty=3)
}

legend("topright",
       legend=c("CR, Epstein-Zin","CR, CRRA",
                "CR, Epstein-Zin & lower damages","CR, CRRA & lower damages",
                "Traeger (2023)"),
       lty=NaN,
       col = c("red","red","red","red",
               "#00000033"),
       pch = c(15,17,0,2,19),
       cex=1,
       lwd=1,
       box.lty=1, box.lwd=1,bg = "white")

dev.off()



data_ACE$r <- data_ACE$prpp + .019
data_ACE$D <- .03*(data_ACE$Damages=="baseline")+.10*(data_ACE$Damages=="HSP")
plot(0,0,col="white",xlim=c(0.015,.045),ylim=c(log(10),log(20000)))
for(i in 1:dim(data_ACE)[1]){
  points(data_ACE$r[i],log(data_ACE$SCC[i]),pch=19,col="#00000033",cex=8*sqrt(data_ACE$D[i]))
}
eq <- lm(SCC~r+I(r^2)+D,data=data_ACE)
eq <- lm(SCC~I(r)+I(D^2),data=data_ACE)
summary(eq)

for(i in 1:length(all_rho)){
  points(all_5yrRates_EZ[i],log(all_scc_EZ[i]),
         pch=19,col="#77000033",cex=8*sqrt(1-model_sol$target_vector["ECumD4"]))
}
for(i in 1:length(all_rho)){
  points(all_5yrRates_EZ_lowDamage[i],log(all_scc_EZ_lowDamage[i]),
         pch=19,col="#77000033",cex=8*sqrt(factor_mult_Damage*(1-model_sol$target_vector["ECumD4"])))
}
for(i in 1:length(all_rho)){
  points(all_5yrRates_CRRA[i],log(all_scc_CRRA[i]),
         pch=19,col="#77000033",cex=8*sqrt(1-model_sol$target_vector["ECumD4"]))
}
for(i in 1:length(all_rho)){
  points(all_5yrRates_CRRA_lowDamage[i],log(all_scc_CRRA_lowDamage[i]),
         pch=19,col="#77000033",cex=8*sqrt(factor_mult_Damage*(1-model_sol$target_vector["ECumD4"])))
}


