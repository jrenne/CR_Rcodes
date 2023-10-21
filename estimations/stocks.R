#Initial model - solves model with initial parameters
source("./procedure/Functions_General_v2.R")

library("viridis") 


# #models for different scenarios
# indic_stocks<-0
# model_new<-model_sol
# model_new$parameters$mu_d<-0
# model_new<-model_solve(model_new,model_new$theta0,
#                        indic_mitig = F,mu.chosen=model_sol$mu)
# model_new2<-model_solve(model_new,model_new$theta0)

indic_stocks<-1
##################
#Data information#
##################
#*X = delc, ytilde, E, Eind, F, Mat, Mup, Mlo, Tat, Tlo, CumD, CumE, Cumdc, H,
#*Cum_div,pd,r.s,Cum_rs,r.1
#*W
##################
#STOCKS&DIVIDENDS#
##################
hs<-model_sol$Tmax
  
model_stock<-stock.solve(model_sol,hs)
EV<-EV.fct(model_stock,hs)

# #test pd
# pd_test<-pd.real.fct(model_stock,hs,5)

