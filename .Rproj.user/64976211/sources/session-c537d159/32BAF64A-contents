# Vector of temperatures for strike K, Digital Option:
K <- c(2,3,3.5)

indic_stocks<-0
#Options on temperatures (atm.)

TAT <- which(model$names.var.X=="T_at")

omega_ZCB <- matrix(0,model_sol$n.X)
omega_T.at <- omega_ZCB
omega_T.at[TAT] <- 1

# Calculations ----
xoptions<-exp(seq(-5,5,length.out = 1000))                                          #grid integral v.hat

q3 <- sequential_hcl(5, "YlOrRd")
q3 <-q3[3:1]

x.lim  <-c(2020,2100)
pos.opt<-c(3,3,3)

P.OpD.names <- NULL
for(k in 1:length(K)){
  eval(parse(text = gsub(" ","",
                         paste("P.OpD.names <- c(P.OpD.names,expression(paste(T[K],' = ',",K[k],",sep='')))")
  )))
}



# P.OpD.k     <-NaN
# P.OpD_RNF.k <-NaN                                                               #Risk-neutrality w/ Fourier

#Normalization with ZCB
Price.ZC <- varphi(model_sol,omega_ZCB,H = H)[[3]]

cl <- makeCluster(number.of.cores)
registerDoParallel(cl)

save.image("outputs/toto.Rdata")
clusterEvalQ(cl,load("outputs/toto.Rdata"))

clusterEvalQ(cl,library(MASS))
clusterEvalQ(cl,library(expm))

all.Probas.Q <- foreach(h = 1:H, .combine=rbind) %dopar% {
  Temperature.Q<-c(varphi.hat.fast(model_sol,omega = omega_ZCB,
                                   H=h,x=xoptions,
                                   a = - omega_T.at,
                                   b = - K)/Price.ZC[h])
  100*Temperature.Q
}

all.Probas.P <- foreach(h = 1:H, .combine=rbind) %dopar% {
  Temperature.P <- 1 - fourier(model_sol,x=xoptions,
                               K,h,TAT)
  100*Temperature.P
}

stopCluster(cl)

k <- 3
plot(all.Probas.P[,k])
lines(all.Probas.Q[,k])

# make sure probas are increasing:
for(k in 1:length(K)){
  cdf.Q <- pmax(pmin(all.Probas.Q[,k],100),0)
  cdf.P <- pmax(pmin(all.Probas.P[,k],100),0)
  for(i in 2:length(cdf.P)){
    if(cdf.Q[i]<cdf.Q[i-1]){
      cdf.Q[i] <- cdf.Q[i-1]
    }
    if(cdf.P[i]<cdf.P[i-1]){
      cdf.P[i] <- cdf.P[i-1]
    }
  }
  all.Probas.Q[,k] <- cdf.Q
  all.Probas.P[,k] <- cdf.P
}

splinefunP1    <-splinefun(model_sol$vec_date[2:(H+1)],all.Probas.P[,1],
                           method="hyman")
splinefunQ1    <-splinefun(model_sol$vec_date[2:(H+1)],all.Probas.Q[,1],
                           method="hyman")
splinefunP2    <-splinefun(model_sol$vec_date[2:(H+1)],all.Probas.P[,2],
                           method="hyman")
splinefunQ2    <-splinefun(model_sol$vec_date[2:(H+1)],all.Probas.Q[,2],
                           method="hyman")
splinefunP3    <-splinefun(model_sol$vec_date[2:(H+1)],all.Probas.P[,3],
                           method="hyman")
splinefunQ3    <-splinefun(model_sol$vec_date[2:(H+1)],all.Probas.Q[,3],
                           method="hyman")
P.OpDPs<-list()
P.OpDQs<-list()
grid   <-model_sol$vec_date[2]:model_sol$vec_date[H+1]

P.OpDPs[[1]]<-splinefunP1(grid)
P.OpDPs[[2]]<-splinefunP2(grid)
P.OpDPs[[3]]<-splinefunP3(grid)

P.OpDQs[[1]]<-splinefunQ1(grid)
P.OpDQs[[2]]<-splinefunQ2(grid)
P.OpDQs[[3]]<-splinefunQ3(grid)


arrow.pos.x0  <-c(model_sol$vec_date[2]+30,
                  model_sol$vec_date[2]+40,
                  model_sol$vec_date[2]+50)
arrow.pos.x1  <-arrow.pos.x0 - 5
arrow.pos.y0<-c(P.OpDQs[[1]][grep(arrow.pos.x0[1],grid)],
                P.OpDQs[[2]][grep(arrow.pos.x0[2],grid)],
                P.OpDQs[[3]][grep(arrow.pos.x0[3],grid)])
arrow.pos.y0 <- (arrow.pos.y0>50)*(arrow.pos.y0-20) +
  (arrow.pos.y0<=50)*(arrow.pos.y0+20)
arrow.pos.y1<-c(P.OpDQs[[1]][grep(arrow.pos.x1[1],grid)],
                P.OpDQs[[2]][grep(arrow.pos.x1[2],grid)],
                P.OpDQs[[3]][grep(arrow.pos.x1[3],grid)])
text.pos.x  <- arrow.pos.x0 + 3
text.pos.y  <- arrow.pos.y0 - 4

# Plot ----
FILE = "/outputs/Figures/Figure_Option_Digital.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=10, width=7, height=4)
par(plt=c(.1,.95,.1,.95))
plot(grid,P.OpDPs[[1]],type="l",
     cex.main=1.0,
     cex.axis=1.0,#cex.lab=1.5,
     col="white",
     lwd=2,lty=2,
     ylab="In percent",xlab="Maturity",
     las=1,
     ylim=c(0,110),xlim=x.lim)
for (k in 1:length(K)){
  lines(grid,P.OpDQs[[k]],
        col=q3[k],
        lwd=3,lty=1)
  lines(grid,P.OpDPs[[k]],
        col=q3[k],
        lwd=3,lty=2)
  arrows(x0=arrow.pos.x0[k],y0=arrow.pos.y0[k],
         x1=arrow.pos.x1[k],y1=arrow.pos.y1[k],
         length = 0.1,col = q3[k],lwd=2)
  text(text.pos.x[k],text.pos.y[k],
       labels = P.OpD.names[k],
       cex = 1.2, pos = pos.opt[k], col = q3[k])
}
abline(h=0,col="grey")
legend("topleft",
       legend=c("Option price","(= risk-adjusted probability)",
                "Option price without risk premium","(= physical probability)"),
       lty=c(1,0,2),
       lwd=rep(2,0,2),
       col=c("black","white","black","white"),
       bty = "n",cex=1.0)
dev.off()
