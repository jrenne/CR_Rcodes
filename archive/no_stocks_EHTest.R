indic_stocks<-0

#EH test
model_EH<-model_sol
model_EH$pi<-rep(list(matrix(0,nrow=length(model_EH$pi[[1]]))),
                 length(model_EH$pi))

model_EH_new<-model_new
model_EH_new$pi<-model_EH$pi

model_EH_new2<-model_new2
model_EH_new2$pi<-model_EH$pi

h_EH<-20

r_EH<-compute_yc(model_EH,h_EH)
r_Q <-compute_yc(model_sol,h_EH)

r_EH_new<-compute_yc(model_EH_new,h_EH,)
r_Q_new <-compute_yc(model_new,h_EH)

r_EH_new2<-compute_yc(model_EH_new2,h_EH)
r_Q_new2 <-compute_yc(model_new2,h_EH)


plot(r_EH[,1],r_EH[,2],type="l",
     ylim=c(range(r_Q[,2],r_EH[,2])))
lines(r_Q[,1],r_Q[,2],col="red")

plot(r_EH_new[,1],r_EH_new[,2],type="l",
     ylim=c(range(r_Q_new[,2],r_EH_new[,2])))
lines(r_Q_new[,1],r_Q_new[,2],col="red")

plot(r_EH_new2[,1],r_EH_new2[,2],type="l",
     ylim=c(range(r_Q_new2[,2],r_EH_new2[,2])))
lines(r_Q_new2[,1],r_Q_new2[,2],col="red")
