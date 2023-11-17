
# Max horizon used in sums:
max.h <- 200

number.of.cores <- 8


all_sigmax   <- c(0,0,model$sigmax,model$sigmax)
all_sigmatau <- c(0,model$sigmatau,0,model$sigmatau)
names.model <-c("Deterministic",
                "Only climate uncertainty",
                "Only productivity uncertainty",
                "Two types of uncertainties")

all_gammas <- c(1.2,2,4)
gamma.names <- paste("gamma = ",all_gammas,sep="")

nb.sigmatau <- length(all_sigmatau)
nb.gammas   <- length(all_gammas)

all_sigmax   <- rep(1,length(all_gammas)) %x% all_sigmax
all_sigmatau <- rep(1,length(all_gammas)) %x% all_sigmatau
all_gammas   <- all_gammas %x% rep(1,nb.sigmatau)

names.model.long <- rep(names.model,nb.gammas)
names.model.with.gammas <- paste(names.model,", gamma = ",all_gammas,sep="")


# Run parallel computations:
cl <- makeCluster(number.of.cores)
registerDoParallel(cl)
save.image("toto.Rdata")
clusterEvalQ(cl,load("toto.Rdata"))
clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))
all_scc <- foreach(i = 1:length(names.model.long), 
                   .combine=rbind) %dopar% {
                     
                     model_aux <- model
                     
                     model_aux$sigmax   <- all_sigmax[i]
                     model_aux$sigmatau <- all_sigmatau[i]
                     model_aux$gamma    <- all_gammas[i]
                     
                     res.SCC <- compute.SCC(model_aux,
                                            X0,Dlambda0,Dc0,h=max.h,indic_decompos = 1)
                     
                     res.SCC$SCC_CO2
                   }

stopCluster(cl)
file.remove("toto.Rdata")

res <- matrix(all_scc,nb.sigmatau,nb.gammas)
rownames(res) <- names.model
colnames(res) <- gamma.names

res <- rbind(res,NaN,
             100*(res[4,]/res[1,]-1))
rownames(res)[dim(res)[1]] <- "Effect of joint uncert., %"

print(res)
