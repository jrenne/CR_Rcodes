
#Information
H       <-model$horiz.2100
nb.simul<-H                                                                 
nb.traj <-10000                                                                

x<-exp(seq(-20,40,length.out = 1000))                                           #grid for proposition 8 (Fourier)
a<-model$omega_T.at                                                             #Options related to T_at in X
b<-2                                                                            #Options pay when a'X < b


#Model-implied EV
EV<-EV.fct(model_sol,H)

#Compare with simulations ----
test <-simul.function(model_sol,nb.simul,nb.traj)
test <-test[-length(test)] #vectors of X removed for this exercise

##Create plots with mean and confidence intervals
#Quantiles
quantiles1<-c(0.025,0.975)
Q<-lapply(test,function(m)
  apply(m,1,function(x)quantile(x,quantiles1,na.rm=TRUE)))
#Variance
sample.var<-lapply(test,function(m)apply(m,1,var))
#Mean of simulations
mean<-lapply(test, function(m)apply(m,1,mean))

sample.upper<-lapply(1:length(mean),function(x)mean[[x]]+2*sample.var[[x]]^(1/2))
sample.lower<-lapply(1:length(mean),function(x)mean[[x]]-2*sample.var[[x]]^(1/2))

bounds.simul<-lapply(1:length(mean),function(x){
  rbind(sample.lower[[x]],sample.upper[[x]])
})
remove(sample.upper, sample.lower)

EX    <-EV$EX[-(model_sol$n.Z+1):-(model_sol$n.Z+model_sol$n.eta)]
bounds<-EV$bounds[-(model_sol$n.Z+1):-(model_sol$n.Z+model_sol$n.eta)]

all.equal(mean,EX)

#Plots with conditional mean and variance
for (i in 1:length(EX)){
  plot(1:H,EX[[i]],type="l",ylab = names(EX[i]), col="red",
       ylim=c(range(bounds.simul[[i]],bounds[[i]],na.rm=TRUE)))
  polygon(c(1:H,H:1),
          c(bounds[[i]][1,],bounds[[i]][2,H:1]),
          col =adjustcolor( "red", alpha.f = 0.15), border = NA)
  lines(1:H,mean[[i]], col="blue")
  polygon(c(1:H,H:1),
          c(bounds.simul[[i]][1,],bounds.simul[[i]][2,H:1]),
          col =adjustcolor( "blue", alpha.f = 0.15), border = NA)
}
#end ----


#Mitigation plot ----
#Approximation of the DICE mitigation process.
mudice<-0
mudice[1:model_sol$Tmax]<-pmin(exp(-1/21*log(0.17)*
                                          (1:(length(model_sol$vec_date)+1)))*
                                            0.17,1)                             
plot(model_sol$mu,type="l",
     ylab="Mitigation rate")
lines(mudice,col="red")
legend("bottomright",legend=c("DICE","Own mitigation"),col=c("red","black"),
       lty=c(1,1))
title(main="Comparison mitigation rate DICE vs our model")
#end ----



stop()


#Test Pricing ----
#Swaps price relative to expected atmospheric temperature path
#*normalized by price of ZCB.
T.swaps<-varphi.tilde(model_sol,model$omega_T.at,H)[[1]]/
   varphi(model_sol,model$omega_ZCB,H)[[3]]
print("Premium between T-swaps and ZCB-swaps")
print(T.swaps-EV$EX$T_at[1:H])
#end ----


#Optimization for nonexplosive model ----
print("Non-Explosive Model post calibration")
non_explosive<-EV.fct(model_sol,200)
plot(non_explosive$date,non_explosive$EX$T_at,type="l",
     ylab="")
title(main="Expected atm. temperature over long horizon")

plot(varphi(model_sol,model$omega_ZCB,200)$r.t,type="l",
     ylab="Interest rate in %")
abline(h=0,col="red")
title(main="Term structure of ZCB i.r. over long horizon")
#end ----

#Damages graphical representation
print("Graphical representation of damages following calibration")
#Damages
xD <- exp(seq(-50,100,length.out = 10000)) 
gammaD<-log(seq(0.2,1,length.out=200))

#Quantiles of interest
all_quantiles <-c(0.1,0.5,0.9)

#Matrix of temperatures in X
T2100<-seq(1,4,by=0.2)
X2100<-matrix(0,model_sol$n.X,length(T2100))
X2100[9,]<-T2100

fCumD<-apply(X2100,2,function(i){
  CumD.Fourier(model_sol,H,i,gammaD,xD)})
print("Fourier Done")
quantCumD<-matrix(0,length(all_quantiles),length(T2100))
for(j in 1:length(all_quantiles)){
  quantCumD[j,]<-apply(fCumD,2,function(i)exp(gammaD[which.min(abs(i-all_quantiles[j]))]))
}

rownames(quantCumD)<-all_quantiles
colnames(quantCumD)<-T2100

plot(T2100,quantCumD["0.5",],type="l",ylim=c(range(quantCumD)))
lines(T2100,quantCumD["0.1",],col="lightpink")
lines(T2100,quantCumD["0.9",],col="lightpink")
segments(1,0.975,4,0.9,col="lightskyblue")  