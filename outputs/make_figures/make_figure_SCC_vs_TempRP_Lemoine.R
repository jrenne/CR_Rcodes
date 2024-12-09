# ==============================================================================
# Figure showing relationship between SCC RP and Temperature RP
# ==============================================================================

# Load baseline model:
aux   <- make.Lemoine.model()
model_baseline <- aux$model
X0    <- aux$X0

H <- 200
nb.replic <- 500

Year <- 80 # 2100 in Year years

percent_chge <- 50

all.damages <- c(
  "G. Lemoine",
  "G. Dell-Jones-Olken",
  "G. Bansal-Kiku-Ochoa",
  "L. Nordhaus-Sztorc",
  "L. Barrage-Nordhaus",
  "L. Howard-Sterner",
  "L. Weitzman",
  "L. Barnett-Brock-Hansen",
  "L. Traeger (ACE-base)",
  "L. Traeger (ACE-HSP)",
  "L. Dietz-Stern"
)
all.colors <- c("#B22222",
                "cornsilk4",
                "#008B8B",
                "#458B00",
                "#CD5555",
                "darkorange2",
                "#838B8B",
                "darkorchid3",
                "darkolivegreen",
                "darkolivegreen",
                "darkgoldenrod3")


all.sigmas <- seq(0,2*model_baseline$sigma,length.out=8)

coef.multip.values <- c(1 - percent_chge/100,1,1 + percent_chge/100)

cl <- makeCluster(number.of.cores)
registerDoParallel(cl)

save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))

clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))

all.res <- foreach(i = 1:length(all.damages), .combine=rbind) %dopar% {
  
  damage.type <- all.damages[[i]]
  
  s.values     <- list((1-percent_chge/100)*model_baseline$s,
                       model_baseline$s,(1+percent_chge/100)*model_baseline$s)
  
  all.SCC   <- NULL
  all.NPV   <- NULL
  all.E.T   <- NULL
  all.E.T.Q <- NULL
  
  for(sigma in all.sigmas){
    
    model_baseline$sigma <- sigma
    
    res.SCC <- compute.SCC.Lemoine(model_baseline,
                                   damage.type = damage.type,
                                   X0,H=H,nb.replic=nb.replic,
                                   s.values = s.values,
                                   coef.multip.values = coef.multip.values,
                                   seed0=123)
    
    all.SCC   <- c(all.SCC,res.SCC$SCC)
    all.NPV   <- c(all.NPV,res.SCC$NPV)
    all.E.T   <- c(all.E.T,res.SCC$E.T[Year])
    all.E.T.Q <- c(all.E.T.Q,res.SCC$E.T.risk.adj[Year])
  }
  
  c((all.SCC - all.NPV)/all.SCC,(all.E.T.Q-all.E.T),all.SCC)
}
stopCluster(cl)
file.remove("outputs/toto.Rdata")

all.SCC.RP <- all.res[,1:length(all.sigmas)]
all.T.RP   <- all.res[,length(all.sigmas)+(1:length(all.sigmas))]




FILE = paste("/outputs/Figures/Figure_SCCvsTRP_Lemoine.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=7, height=5.5)

load(file="outputs/results/SCC_vs_TRP_CR.Rdat") # Load results from CR model

par(mfrow=c(1,1))
par(plt=c(.15,.95,.15,.95))
plot(0,0,
     ylim=c(min(-.6,min(all.SCC.RP,all.SCC.RP.CR)),
            max(+.3,max(all.SCC.RP,all.SCC.RP.CR))),
     xlim=c(min(-.1,min(all.T.RP,all.T.RP.CR)),
            max(+.1,max(all.T.RP,all.T.RP.CR))),
     # xlim=c(-.3,1),
     # ylim=c(-10,1),
     col="white",
     xlab="2100 Temperature risk premium, in °C",
     ylab="SCC minus NPV of benefits (as a fraction of SCC)")
abline(h=0,col="grey",lty=3)
abline(v=0,col="grey",lty=3)

for(i in 1:length(all.damages)){
  lines(all.T.RP[i,],all.SCC.RP[i,],col=all.colors[i],lwd=2)
  points(all.T.RP[i,1],all.SCC.RP[i,1],col=all.colors[i],lwd=2,pch=i-1)
}

for(i in 1:dim(all.T.RP.CR)[1]){
  lines(all.T.RP.CR[i,],all.SCC.RP.CR[i,],lty=i+1,lwd=2)
  points(all.T.RP.CR[i,1],all.SCC.RP.CR[i,1],lwd=2,pch=0)
}

legend("topleft",
       legend=c(extract(all.damages,1),
                expression(paste("CR, EZ, ",gamma," = 7",sep="")),
                expression(paste("CR, EZ, ",gamma," = 1.001",sep="")),
                expression(paste("CR, CRRA, ",gamma," = 1.45",sep=""))),
       lty=c(rep(1,length(all.damages)),2:4),
       cex=1,
       col=c(all.colors,rep("black",3)),
       pch=c(0:(length(all.damages)-1),rep(NaN,dim(all.T.RP.CR)[1])),
       lwd=2,bty = "n")

dev.off()
