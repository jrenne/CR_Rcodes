# ==============================================================================
# Figure showing sensitivity of SCC to uncertainty
# ==============================================================================

#Climate beta
y.lim   <- c(-.3,1.2)
y.limscc<- c(0,1000)

# Considered values of mu_N (select three values):
values.of.mu_N <- c(0.001,round(model_sol$parameters$mu_N),
                    1.5*round(model_sol$parameters$mu_N))

# Considered values of mu_D:
seq.mu_D<-seq(0.001,
              4.5*model_sol$parameters$mu_D,
              #0.15,
              len=16)

# Run parallel computations:
cl <- makeCluster(number.of.cores)
registerDoParallel(cl)
save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))
clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))
all_scc <- foreach(i = 1:length(seq.mu_D), 
                   .combine=rbind) %dopar% {
                     
                     model_new <- model_sol
                     model_new$parameters$mu_D <- seq.mu_D[i]
                     model_sol_new <- model_solve(model_new)
                     scc <- scc.fct(model_sol_new,0)
                     
                     scc
                   }

stopCluster(cl)
file.remove("outputs/toto.Rdata")



FILE = paste("/outputs/Figures/Figure_SCC.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=9,width=5, height=3)
par(mfrow=c(1,1))
par(plt=c(.15,.95,.2,.95))
plot(seq.mu_D,all_scc,type="l",lwd=2,col="black",
     ylim=y.limscc,
     xlab=expression(paste("Damage uncertainty, ",mu[D],sep="")),
     main="",
     ylab=expression(paste("SCC (in U.S. $ per ton of ",CO[2],")"),sep=""),
     las=0)
abline(h=0,lty=3,col="black")
abline(v=model_sol$parameters$mu_D,lty=3,col="lightslategrey")
dev.off()

