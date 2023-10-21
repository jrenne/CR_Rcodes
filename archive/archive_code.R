#######################################
#####Check Model2Param/Param2Model#####
#######################################
#mp<-matrix(unlist(Param2Model(Model2Param(model),model)$parameters),ncol=1)
#print(max(abs(mp-matrix(unlist(param),ncol=1))))


######################################
#####Test simulations for Cum_div#####
######################################

# test <-simul.function(model_sol_stock,model_sol_stock$Tmax,nb.traj)
# test <-test[-length(test)]
# ##Create plots with mean and confidence intervals
# #Quantiles
# quantiles1<-c(0.025,0.975)
# Q<-lapply(test,function(m)
#   apply(m,1,function(x)quantile(x,quantiles1,na.rm=TRUE)))
# #Variance
# sample.var<-lapply(test,function(m)apply(m,1,var))
# #Mean of simulations
# mean<-lapply(test, function(m)apply(m,1,mean))[1:15]
# 
# sample.upper<-lapply(1:length(mean),function(x)mean[[x]]+2*sample.var[[x]]^(1/2))
# sample.lower<-lapply(1:length(mean),function(x)mean[[x]]-2*sample.var[[x]]^(1/2))
# 
# bounds.simul<-lapply(1:length(mean),function(x){
#   rbind(sample.lower[[x]],sample.upper[[x]])
# })
# remove(sample.upper, sample.lower)
# 
# EX    <-EV$EX[-(model_sol$n.Z+1):-length(EV$EX)]
# bounds<-EV$bounds[-(model_sol$n.Z+1):-length(EV$bounds)]
# #Plot Cum_div
# plot(1:length(model_sol$vec_date),EX[[15]],type="l",ylab = names(EX[15]),
#      ylim=c(range(bounds.simul[[15]],bounds[[15]],na.rm=TRUE)))
# polygon(c(1:length(model_sol$vec_date),length(model_sol$vec_date):1),
#         c(bounds[[15]][1,],bounds[[15]][2,length(model_sol$vec_date):1]),
#         col =adjustcolor( "red", alpha.f = 0.15), border = NA)
# lines(mean[[15]], col="blue")
# lines(EX[[13]], col="blue")
# polygon(c(1:length(model_sol$vec_date),length(model_sol$vec_date):1),
#         c(bounds.simul[[15]][1,],bounds.simul[[15]][2,length(model_sol$vec_date):1]),
#         col =adjustcolor( "blue", alpha.f = 0.15), border = NA)