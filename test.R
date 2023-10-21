library(viridis)

optmitig.fct<-function(model_sol,nb.simul,nb.traj,theta){
  X.simul<-simul.function(model_sol,nb.simul,nb.traj)$X[[nb.simul]]
  mitig.opt<-apply(X.simul,2,function(x){
    unlist(res.optim(model_sol,theta,Tend=model_sol$Tmax-nb.simul,X=x))})
  
  avg <-apply(mitig.opt,1,mean)
  mu.t<-c(model_sol$mu[1:nb.simul],
          pmin(exp(log(model_sol$mu[nb.simul])+
                     abs(avg[2])*(1:(model_sol$Tmax-nb.simul))),1))
  print(nb.simul)
  return(list(mu.t,avg))
}

# f<-lapply(sample,function(x)c(model_sol$mu[1:x],
#   pmin(-abs(testmit[[x]][[2]][1])+abs(testmit[[x]][[2]][2])*(1:(model_sol$Tmax-x)),1)))
# 
# plot(model_sol$vec_date,head(model_sol$mu,-1),
#      type="l",ylim=c(range(f)))
# for(i in sample){
#   lines(model_sol$vec_date,
#         head(f[[i]],-1),
#         col=col[i])
#   #locator(1)
# }

sample  <-1:50
testmit <-lapply(sample,function(i){optmitig.fct(model_sol,i,50,theta0[[1]])})
col     <-plasma(length(sample))

FILE = paste("/outputs/Figures/Figure_SimulMitig.pdf",sep="")
pdf(file=paste(getwd(),FILE,sep=""),pointsize=7,width=8, height=6)
par(mfrow=c(1,1))
plot(model_sol$vec_date,head(model_sol$mu,-1),
     type="l",ylim=c(range(0,1)))
for(i in sample){
  lines(model_sol$vec_date,
        head(testmit[[i]][[1]],-1),
        col=col[i])
  }
dev.off()

plot(model_sol$vec_date,head(model_sol$mu,-1),
     type="l",ylim=c(range(0,1)))
for(i in sample){
  lines(model_sol$vec_date,
        head(testmit[[i]][[1]],-1),
        col=col[i])
  locator(1)
}