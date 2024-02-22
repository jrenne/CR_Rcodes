

H <- 200
nb.replic <- 500

Year <- 80 # 2100 in Year years

percent_chge <- 50

damage.type <- "G. Lemoine"

all.damage <- list()

all.damage[[1]] <- list(
  damage.type = "G. Lemoine",
  alpha = .0018,
  coef.multip = 1)

all.damage[[2]] <- list(
  damage.type = "G. Dell-Jones-Olken",
  alpha = .00137,
  coef.multip = 1)

all.damage[[3]] <- list(
  damage.type = "G. Nordhaus-Sztorc",
  alpha = .00026,
  coef.multip = 1)

all.damage[[4]] <- list(
  damage.type = "G. Weitzman",
  alpha = c(.000075,3.25),
  coef.multip = 1)

all.damage[[5]] <- list(
  damage.type = "L. Nordhaus-Sztorc",
  alpha = .00266,
  coef.multip = 1)

all.damage[[6]] <- list(
  damage.type = "L. Weitzman",
  alpha = c(20.64,6.081,6.754),
  coef.multip = 1)

all.damage[[7]] <- list(
  damage.type = "L. Barnett-Brock-Hansen",
  alpha = c(0.0127,0.0005,.0005),
  coef.multip = 1)

# all.damage[[8]] <- list(
#   damage.type = "L-Bansal",
#   alpha = c(.0050,.0033,.0011,.0011),
#   coef.multip = 1)

# Baseline model (Lemoine, 2021):
mu0    <- .047
mu1    <- .0134
sigma  <- .014
l0     <- .022
l1     <- .012
psi1   <- .2
psi2   <- .393
psi    <- .0023
nu     <- 5.35
Mpre   <- 605
gamma0 <- .28
gamma1 <- .015
phi    <- .0394
S      <- 3.13
s      <- nu * log(2) / S

eta    <- 1.45
beta   <- 1 - .015

all.sigmas <- seq(0,2*sigma,length.out=8)

model <- list(
  eta = eta,
  beta = beta,
  mu0 = mu0,
  mu1 = mu1,
  sigma = sigma,
  l0 = l0,
  l1 = l1,
  psi1 = psi1,
  psi2 = psi2,
  psi = psi,
  nu = nu,
  Mpre = Mpre,
  gamma0 = gamma0,
  gamma1 = gamma1,
  phi = phi,
  s = s,
  damage = all.damage[[which(extract(all.damage,1)==damage.type)]]
)

# Initial values:
C0  <- 55 * 1.1772 # in trillion USD 2014
L0  <- 7.1 * 10^9
T0  <- .8794
M10 <- 706
M20 <- 148
X0 <- c(C0,L0,T0,M10,M20)


coef.multip.values <- c(1 - percent_chge/100,1,1 + percent_chge/100)

cl <- makeCluster(number.of.cores)
registerDoParallel(cl)

save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))

clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))

all.res <- foreach(i = 1:length(all.damage), .combine=rbind) %dopar% {
  
  model$damage <- all.damage[[i]]
  
  s.values     <- list((1-percent_chge/100)*model$s,
                       model$s,(1+percent_chge/100)*model$s)
  
  all.SCC   <- NULL
  all.NPV   <- NULL
  all.E.T   <- NULL
  all.E.T.Q <- NULL
  
  for(sigma in all.sigmas){
    
    model$sigma <- sigma
    
    res.SCC <- compute.SCC.Lemoine(model,X0,H=H,nb.replic=nb.replic,
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
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10,width=7, height=4)

par(mfrow=c(1,1))
par(plt=c(.15,.95,.2,.95))
plot(0,0,
     ylim=c(-.7,max(.3,max(all.SCC.RP))),
     xlim=c(-.2,max(.2,max(all.T.RP))),
     col="white",
     xlab="2100 Temperature risk premium, in °C",
     ylab="SCC minus NPV of benefits (as a fraction of SCC)")
abline(h=0,col="grey",lty=3)
abline(v=0,col="grey",lty=3)

for(i in 1:length(all.damage)){
  lines(all.T.RP[i,],all.SCC.RP[i,],col=i+1,lwd=2)
  points(all.T.RP[i,1],all.SCC.RP[i,1],col=i+1,lwd=2,pch=0)
}

load(file="outputs/results/SCC_vs_TRP_CR.Rdat")
for(i in 1:dim(all.T.RP.CR)[1]){
  lines(all.T.RP.CR[i,],all.SCC.RP.CR[i,],lty=i,lwd=2)
  points(all.T.RP.CR[i,1],all.SCC.RP.CR[i,1],lwd=2,pch=0)
}

legend("topleft",
       legend=c(extract(all.damage,1),
                expression(paste("CR, EZ, ",gamma," = 7",sep="")),
                expression(paste("CR, EZ, ",gamma," = 1.001",sep="")),
                expression(paste("CR, CRRA, ",gamma," = 1.45",sep=""))),
       lty=c(rep(1,length(all.damage)),1:3),
       cex=1,
       col=c(1+(1:length(all.damage)),rep("black",3)),
       lwd=2,bty = "n")

dev.off()
