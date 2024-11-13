
#clear environment
rm(list=ls(all=T)) 

seed <- rnorm(1)

nb_discret_values <- 10

source("set.of.proc.Lemoine.R")

indic.uncertainty.alpha <- 1
indic.uncertainty.S     <- 1
indic.uncertainty.C     <- 1

H <- 200
nb.replic <- 500

damage.type <- "Lemoine"
# damage.type <- "G-N"
# damage.type <- "G-W"
# damage.type <- "L-N"
# damage.type <- "L-W"
# damage.type <- "G-DJO"
# damage.type <- "L-Barnett"
# damage.type <- "L-Bansal"

mu0    <- .047
mu1    <- .0134
alpha  <- .0018
sigma  <- .014 * indic.uncertainty.C
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

#M0     <- 800 # where radiative forcings are linearized

s <- nu * log(2) / S

C0  <- 55 * 1.1772 # in trillion USD 2014
L0  <- 7.1 * 10^9
T0  <- .8794
M10 <- 706
M20 <- 148
# NB: M = M10 + M20, and divide by 2.13 to have M expressed in ppmv CO2

eta <- 1.45
beta <- 1 - .015

X0 <- c(C0,L0,T0,M10,M20)

model <- list(
  mu0 = mu0,
  mu1 = mu1,
  alpha = alpha,
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
  #M0 = M0,
  s = s
)

#nb.replic <- 4
seed <- seed + 1
set.seed(seed)
res.simul     <- simul.model.Lemoine(model,X0,H=H,nb.replic = nb.replic,
                             damage.type = damage.type)


res.mean.traj <- lapply(res.simul,function(x){apply(x,1,mean)})

par(mfrow=c(2,3))
plot(res.simul$all.C[,1],type="l")
lines(res.mean.traj$all.C,col="blue",lty=2)
plot(res.simul$all.C[,1]/res.simul$all.L[,1]*10^9,type="l")
lines(res.mean.traj$all.C/res.mean.traj$all.L*10^9,col="blue",lty=2)
plot(res.simul$all.T[,1],type="l")
lines(res.mean.traj$all.T,col="blue",lty=2)
plot(res.simul$all.M1[,1],type="l",
     ylim=c(0,max(res.simul$all.M1,res.simul$all.M2)))
lines(res.mean.traj$all.M1,col="blue",lty=2)
lines(res.simul$all.M2[,1],type="l",col="red")
lines(res.mean.traj$all.M2,col="blue",lty=2)


#stop()


# SCC

sumE <- function(x,y=matrix(1,dim(x)[1],1)){
  return(sum(apply(x,1,mean)*apply(y,1,mean)))
}


matrix.beta.h <- matrix(exp(-(1-beta)*(1:H)),H,nb.replic)

seed <- seed + 1

set.seed(seed)
res.simul0 <- simul.model.Lemoine(model,X0,H=H,nb.replic = nb.replic,
                          damage.type = damage.type)
C.L0 <- res.simul0$all.C/res.simul0$all.L
U   <- matrix.beta.h * C.L0^(1-eta)/(1-eta) #* res.simul1$all.L/L0 # difference with Lemoine
U0  <- sumE(U)

set.seed(seed)
X1 <- X0
shock <- 1
X1[4] <- X0[4] - shock * model$psi1
X1[5] <- X0[5] - shock * (1 - model$psi1) * model$psi2

res.simul1 <- simul.model.Lemoine(model,X0=X1,H=H,nb.replic = nb.replic,
                          damage.type = damage.type)
C.L1 <- res.simul1$all.C/res.simul1$all.L
U   <- matrix.beta.h * C.L1^(1-eta) /(1-eta) #* res.simul1$all.L/L0 # difference with Lemoine
U1  <- sumE(U)

SCC <- (U1 - U0)/((C0/L0)^(-eta)/L0) * 10^12/10^9 / 3.667
print(SCC)

M0 <- matrix.beta.h * (C.L0)^(-eta)/(C0/L0)^(-eta)
M1 <- matrix.beta.h * (C.L1)^(-eta)/(C0/L0)^(-eta)
B <- C.L1 - C.L0
SCC.check <- L0 * sumE(M0*B) * 10^12/10^9 / 3.667
SCC.check <- L0 * sumE(M1*B) * 10^12/10^9 / 3.667
MM <- .5*(M0+M1)
SCC.check <- L0 * sumE(MM*B) * 10^12/10^9 / 3.667
print(SCC.check)

NPV <- L0 * sumE(MM,B) * 10^12/10^9 / 3.667

cbind(SCC,SCC.check,NPV)

# # Change in wealth:
# L0 * sum(apply(M0*C.L1-M0*C.L0,1,mean)) * 10^12/10^9 / 3.667
# L0 * sum(apply(M1*C.L1-M1*C.L0,1,mean)) * 10^12/10^9 / 3.667
# L0 * sum(apply(M1*C.L1-M0*C.L0,1,mean)) * 10^12/10^9 / 3.667


# Lemoine definition of SCC (eq.5):
epsilon <- (res.simul1$all.C - res.simul0$all.C)/res.simul0$all.C
SCC.Lemoine <- 
  sumE(matrix.beta.h * (C.L0/(C0/L0))^(-eta) * epsilon * res.simul0$all.C) * 
  10^12/10^9 / 3.667
print(SCC.Lemoine)


#stop()



# Climate sensitivity:
x <- seq(0,10,by=.01)
plot(x,dlnorm(x,meanlog = 1.10704,sdlog = 0.264),type="l")
S.values <- qlnorm(seq(.05,.95,length.out=nb_discret_values),meanlog = 1.10704,sdlog = 0.264)
# Values of s:
s.values <- nu*log(2)/S.values

# alpha:
x <- seq(0,.0094,length.out=500)
#x <- seq(0,.00094,length.out=500)
alpha_max <- .5/exp(4.2882+0.1494^2/2)
x <- seq(0,alpha_max,length.out=500)
plot(x,dlnorm(x,meanlog = -6.7342,sdlog = 1.4684),type="l")
xxx <- plnorm(alpha_max,meanlog = -6.7342,sdlog = 1.4684)
xx <- .01
alpha.values <- qlnorm(seq(xx,xxx-xx,length.out=nb_discret_values),
                       meanlog = -6.7342,sdlog = 1.4684)
#alpha.values <- qlnorm(seq(.3,.7,length.out=5),meanlog = -6.7342,sdlog = 1.4684)

# Damages:
x <- seq(0,1,length.out=500)
plot(x,dlnorm(x,meanlog = -2.446,sdlog = 1.476),type="l")

# Sum T until + 50 yrs
x <- seq(0,200,length.out=500)
plot(x,dlnorm(x,meanlog = 4.2882,sdlog = .1494),type="l")

sumT <- rlnorm(10000,meanlog = 4.2882,sdlog = .1494)
Alpha <- rlnorm(10000,meanlog = -6.7342,sdlog = 1.4684)

Theta <- sumT*alpha
#mean(log(Theta))
#sd(log(Theta))
#mean(sumT)

# Values of alpha:
#alpha.values <- seq(.001,0.0094,length.out=5)

seed0 <- 123
seed  <- seed0  

all.sim.C0 <- NULL
all.sim.C1 <- NULL
all.sim.L  <- NULL
all.sim.T  <- NULL
all.SCC <- NULL

if(indic.uncertainty.S==0){
  s.values <- s
}
if(indic.uncertainty.alpha==0){
  alpha.values <- alpha
}


for(s in s.values){
  for(alpha in alpha.values){
    
    model <- list(
      mu0 = mu0,
      mu1 = mu1,
      alpha = alpha,
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
      M0 = M0,
      s = s
    )
    
    set.seed(seed0)
    res.simul0 <- simul.model.Lemoine(model,X0,H=H,nb.replic = nb.replic,
                              damage.type = damage.type)
    set.seed(seed0)
    X1 <- X0
    X1[4] <- X0[4] - shock * model$psi1
    X1[5] <- X0[5] - shock * (1 - model$psi1) * model$psi2
    res.simul1 <- simul.model.Lemoine(model,X0=X1,H=H,nb.replic = nb.replic,
                              damage.type = damage.type)
    
    matrix.beta.h <- matrix(exp(-(1-beta)*(1:H)),H,nb.replic)
    
    C.L0 <- res.simul0$all.C/res.simul0$all.L
    U    <- matrix.beta.h * C.L0^(1-eta)/(1-eta)
    U0   <- sumE(U)
    
    C.L1 <- res.simul1$all.C/res.simul0$all.L
    U   <- matrix.beta.h * C.L1^(1-eta)/(1-eta)
    U1  <- sumE(U)
    SCC <- (U1 - U0)/((C0/L0)^(-eta)/L0) * 10^12/10^9 / 3.667
    all.SCC <- c(all.SCC,SCC)
    
    all.sim.C0 <- cbind(all.sim.C0,res.simul0$all.C)
    all.sim.C1 <- cbind(all.sim.C1,res.simul1$all.C)
    all.sim.L  <- cbind(all.sim.L,res.simul0$all.L)
    all.sim.T  <- cbind(all.sim.T,res.simul0$all.T)
  }
}

matrix.beta.h <- matrix(exp(-(1-beta)*(1:H)),H,dim(all.sim.C0)[2])

C.L0 <- all.sim.C0/all.sim.L
U    <- matrix.beta.h * C.L0^(1-eta)/(1-eta)
U0   <- sumE(U)

C.L1 <- all.sim.C1/all.sim.L
U   <- matrix.beta.h * C.L1^(1-eta)/(1-eta)
U1  <- sumE(U)

SCC <- (U1 - U0)/((C0/L0)^(-eta)/L0) * 10^12/10^9 / 3.667
print(SCC)

B <- C.L1 - C.L0
M0 <- matrix.beta.h * (C.L0)^(-eta)/(C0/L0)^(-eta)
M1 <- matrix.beta.h * (C.L1)^(-eta)/(C0/L0)^(-eta)
MM <- .5*(M0+M1)
SCC.check <- L0 * sumE(MM*B) * 10^12/10^9 / 3.667

NPV <- L0 * sumE(MM,B) * 10^12/10^9 / 3.667
cbind(SCC,SCC.check,NPV)

# Compute changes in wealth:
chge.wealth <- L0 * sumE(M1*C.L1-M0*C.L0) * 10^12/10^9 / 3.667

# Apply Lemoine (2021) decomposition:
epsilon <- (all.sim.C1 - all.sim.C0)/all.sim.C0

SCC.Lemoine <- 
  sumE(matrix.beta.h * (C.L0/(C0/L0))^(-eta) * epsilon * all.sim.C0) * 
  10^12/10^9 / 3.667

part.determ <- sum(
  matrix.beta.h[,1] / C0^(-eta) * (all.sim.L[,1]/L0)^eta *
    apply(all.sim.C0,1,mean)^(1-eta) * apply(epsilon,1,mean) * 
    10^12/10^9 / 3.667)

part.precaut <- sum(
  .5 * eta * (eta + 1) *
    matrix.beta.h[,1] / C0^(-eta) * (all.sim.L[,1]/L0)^eta *
    apply(all.sim.C0,1,mean)^(-1-eta) * apply(epsilon,1,mean) * 
    apply(all.sim.C0,1,var) *
    10^12/10^9 / 3.667)

part.Dscaling <- sum(
  - eta * 
    matrix.beta.h[,1] / C0^(-eta) * (all.sim.L[,1]/L0)^eta *
    apply(all.sim.C0,1,mean)^(-1-eta) * apply(epsilon,1,mean) * 
    apply(all.sim.C0,1,var) *
    10^12/10^9 / 3.667)

part.insur <- sum(
    matrix.beta.h[,1] / C0^(-eta) * (all.sim.L[,1]/L0)^eta *
    (apply(all.sim.C0^(1-eta)*epsilon,1,mean) - 
       apply(all.sim.C0^(1-eta),1,mean)*apply(epsilon,1,mean)) *
    10^12/10^9 / 3.667)


part.123 <- sum(
  matrix.beta.h[,1] / C0^(-eta) * (all.sim.L[,1]/L0)^eta *
    apply(all.sim.C0^(1-eta),1,mean) * apply(epsilon,1,mean) * 
    10^12/10^9 / 3.667)

SCC.check <- part.determ + part.precaut + part.Dscaling + part.insur
  
print(cbind(
  SCC.Lemoine,part.determ,part.precaut,
  part.Dscaling,part.insur,SCC.check,part.123)
)

part.123.check <- part.determ + part.precaut + part.Dscaling
print(cbind(part.123,part.123.check))


# Compute Temperature risk premium:
T <- all.sim.T
T.risk.adj <- apply(MM*T,1,mean)/apply(MM,1,mean)
plot(T.risk.adj,type="l")
lines(apply(T,1,mean),col="red")


