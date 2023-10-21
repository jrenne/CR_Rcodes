load(file="estimations/solve_param_estim/param_final_est.Rdat")

#model_est$parameters$mu_d    <- .03
#model_est$parameters$mu_n    <- .000001
#model_est$parameters$gamma   <- .5
# model_est$parameters$sigma_a <- 0.0000001
# model_est$bounds$min_sigma_a <- .00000000001
# model_est$parameters$delta <- .9
#model_est$parameters$b_sk <- 0
#model_est$parameters$delta <- (1 - .015)^5

#model_sol$X[13]

#RESET FILTER
FILTER                        <- rep(0,length(FILTER))
FILTER[1]                     <- 1 # A_bar
FILTER[2]                     <- 1 # sigma_a
FILTER[6+n.eta]               <- 1 # ell.0.N
FILTER[15+n.eta]              <- 1 # ell.1.N[9]
FILTER[7+n.eta+n.Z+n.W]       <- 1 # rho.N
FILTER[17+n.eta+n.Z+n.W]      <- 1 # ell.1.D[9]
FILTER[9+n.eta+2*(n.Z+n.W)]   <- 1 # mu.D
FILTER[10+n.eta+2*(n.Z+n.W)]  <- 1 # mu.N
FILTER[38+n.eta+2*(n.Z+n.W)]  <- 1 # pback
FILTER[40+n.eta+2*(n.Z+n.W)]  <- 1 # sigma_eta_f
FILTER[43+n.eta+2*(n.Z+n.W)]  <- 1 # a_sat
FILTER[68+n.eta+3*(n.Z+n.W)]  <- 1 # c_sat

param_vector <- Model2Param(model_est)
param_est    <- param_vector[FILTER==1]

mim<-compute.moments(param_est,model_est,FILTER,horiz,theta0)
print("New moments and targets achieved")
print(cbind(target_vector,mim))
print(cbind(c("mu_n","mu_d","rho.N","ell0.N","ell1.N","ell1.D","pback",
              "a_sat","c_sat"),
            c(sprintf(format.nb,model_est$parameters$mu_n),
              sprintf(format.nb,model_est$parameters$mu_d),
              sprintf(format.nb,model_est$parameters$rho.N),
              sprintf(format.nb,model_est$parameters$ell0.N),
              sprintf(format.nb,model_est$parameters$ell1.N[9]),
              sprintf(format.nb,model_est$parameters$ell1.D[9]),
              sprintf(format.nb,model_est$parameters$pback),
              sprintf(format.nb,model_est$parameters$a_sat),
              sprintf(format.nb,model_est$parameters$c_sat)
            )))


model_sol<-model_solve(model_est,theta0,indic_mitig = TRUE)

mu_u <- mu_u.t.fct.all(model_sol)
h <- 1
mu_u$mu_u1.t1[[h]][6]

h1 <- 1
h2 <- 10
cbind(mu_u$mu_u1.t1[[h1]],mu_u$mu_u1.t1[[h2]])



model_sol$X[13] <- 0

H <- 200

alpha_c <- model_sol$omega_ZCB
alpha_c[13] <- 1

Delta <- model_sol$omega_ZCB
Delta[6] <- -1

# Delta <- model_sol$omega_ZCB
# Delta[5] <- -1

# Delta <- model_sol$omega_ZCB
# Delta[9] <- -1

model_shocked <- model_sol
model_shocked$X <- model_shocked$X + Delta

baseline <- varphi(model_sol,alpha_c,H=H)
shocked  <- varphi(model_shocked,alpha_c,H=H)

# check utility:
sum(baseline$P.t)*model_sol$c0
exp(model_sol$u0)/model_sol$c0


P.baseline <- sum(baseline$P.t)
P.shocked  <- sum(shocked$P.t)

print(P.shocked - P.baseline)


cbind(baseline$P.t,shocked$P.t)


ZC <- varphi(model_sol,model_sol$omega_ZCB,H=H)
plot(ZC$P.t,type="l")
plot(ZC$r.t,type="l")




#stop()



EV <- EV.fct(model_sol,h=H)
EV_shocked <- EV.fct(model_shocked,h=H)

cbind(EV$EX$Cum_dc,EV_shocked$EX$Cum_dc,EV_shocked$EX$Cum_dc-EV$EX$Cum_dc)

conso.chge <- EV_shocked$EX$Cum_dc-EV$EX$Cum_dc

all.EC.baseline <- NULL
all.EC.shocked  <- NULL
U    <-matrix(0,model_sol$n.X,1)
U[13,] <- 1

all.P.0 <- NULL
all.P.1 <- NULL
all.P.2 <- NULL
all.PPP <- NULL

for(h in 1:H){
  E_C_h_baseline <- psi(model_sol,u=1,h,i=13,X=model_sol$X)
  E_C_h_shocked  <- psi(model_sol,u=1,h,i=13,X=model_shocked$X)
  
  res.baseline <- multi.lt.fct(model_sol,U,h,X=model_sol$X,t=0)
  res.shocked  <- multi.lt.fct(model_sol,U,h,X=model_shocked$X,t=0)
  
  all.EC.baseline <- c(all.EC.baseline,E_C_h_baseline)
  all.EC.shocked  <- c(all.EC.shocked,E_C_h_shocked)
  
  P.0 <-  model_sol$c0 * ZC$P.t[h] * res.baseline$uX_t.h * (c(Delta) %*% res.baseline$psi1)
  P.1 <-  model_sol$c0 * ZC$P.t[h] * res.baseline$uX_t.h * (c(Delta) %*% ZC$varphi1[[h]])
  P.2 <-  model_sol$c0 * baseline$P.t[h] * (c(Delta) %*% baseline$varphi1[[h]]) -
    P.0 - P.1
  
  all.PPP <- c(all.PPP,
               model_sol$c0 * baseline$P.t[h] * (c(Delta) %*% baseline$varphi1[[h]]))
  
  all.P.0 <- c(all.P.0,P.0)
  all.P.1 <- c(all.P.1,P.1)
  all.P.2 <- c(all.P.2,P.2)
}

plot(all.P.0,type="l",ylim=c(min(all.P.0,all.P.1,all.P.2),
                             max(all.P.0,all.P.1,all.P.2)))
lines(all.P.1,col="blue")
lines(all.P.2,col="red")

sum(all.P.0 + all.P.1 + all.P.2)

sum(all.P.0 + all.P.2)


conso.chge <- all.EC.shocked - all.EC.baseline

r.t.Z.spread <- ZC$r.t + 1.67
P.t.Z.spread <- exp(-r.t.Z.spread/100*(1:H)*model_sol$tstep)

NPV.conso.chge          <- sum(conso.chge*ZC$P.t) * model_sol$c0 * 10^12/10^9
NPV.conso.chge.Z.spread <- sum(conso.chge*P.t.Z.spread) * model_sol$c0 * 10^12/10^9

plot(EV$EX$Cum_dc,type="l")
lines(EV_shocked$EX$Cum_dc,col="red")

#cbind(model_sol$X,model_shocked$X)

scc <- - model_sol$c0 * mu_u.t.fct(model_sol)[[2]][6] * 10^(3)/
  (1-model_sol$parameters$delta)

print(scc)
print(NPV.conso.chge)
print(NPV.conso.chge/scc)
print(NPV.conso.chge.Z.spread)


HH <- 95

FILE = paste("/outputs/Figures/Figure_SCC_Zspread.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=6, height=5)

par(mfrow=c(2,1))
par(plt=c(.15,.95,.2,.8))
plot(model_sol$tstep*(1:HH),conso.chge[1:HH] * model_sol$c0 * 10^12/10^9,
     type="l",lwd=2,ylab="Expected benefit (in dollars)",
     xlab="horizon, in years",las=1,
     main="(a) Expected benefits (- 1 ton of carbon)")

plot(model_sol$tstep*(1:HH),ZC$r.t[1:HH],lwd=2,type="l",
     ylab="interest rates (in percent)",
     xlab="maturity, in years",
     las=1,
     main="(b) Discount rates",
     ylim=c(0,max(r.t.Z.spread)))
lines(model_sol$tstep*(1:HH),r.t.Z.spread[1:HH],lwd=3,lty=3)
dev.off()


stop()

HH <- 100

# for(i in 1:98){
#   model_sol$pi[[i]] <- 0*model_sol$pi[[i]]
# }

res1 <- varphi(model_sol,alpha_c,H=HH)
res0 <- multi.lt.fct(model_sol,alpha_c,h=HH)
ZC   <- varphi(model_sol,model_sol$omega_ZCB,H=HH)

# Price of C_{t+h}/E_t(C_{t+h}):
p0 <- res1$P.t[HH]/res0$uX_t.h
# Price of ZC (same expected payoff)
p1 <- ZC$P.t[HH]
cbind(p0,p1,p1/p0)


stop()


# Infini:
print(model_sol$mu_u0 + t(model_sol$mu_u1)%*%model_sol$X)

# Optimized utility at 0:
print(model_sol$u0 - log(model_sol$c0))
print(exp(model_sol$u0 - log(model_sol$c0)))

# reproduce optimized utility at 0:
res <- mu_u.t.fct(model_sol)
print(res$mu_u0.t + t(res$mu_u1.t)%*%model_sol$X)

P.mu_u.t1 <-mu_u.t.fct.all(model_sol)
i <- 1
print(P.mu_u.t1$mu_u0.t1[i] + t(unlist(P.mu_u.t1$mu_u1.t1[i]))%*%model_sol$X)



# Check mu:
opt <- res.optim(model_sol,c(-10,0))
theta <- opt[[1]]
plot(mu.function(model_sol,theta))

delta <- model_sol$parameters$delta
dc <- model_sol$omega0.inf[1]
M <- delta * exp(-dc)
M^2

baseline$P.t[1]
baseline$P.t[2]

# Pbm with the price of C_{t+h}, which should be equal to delta^h.

