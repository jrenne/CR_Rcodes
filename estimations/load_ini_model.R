# ==============================================================================
# Initial model - solves model with initial parameters
# ==============================================================================

#number of periods before stationary model:
Tmax   <- 100

#number of years in each period t \in(1:Tmax):
tstep  <- 5

#vector of dates
vec_date <- seq(2020,by=tstep,len=Tmax-1)

#Initial values of theta_a/_b for mitig optim. ~ DICE2016
#theta0  <- c(log(0.17),-1/21*log(0.17))
theta0  <- c(-2,-.1)

#Determine of max of iterations in mitig optim.
MAXIT   <- 200

#Define the horizon for optimization
horiz   <- (2100-vec_date[1])/tstep


#---- state vector -------------------------------------------------------------

#---- INITIALIZATION ----
#The World in Data + DICE
eps_0 <- 5.9  # DICE2023
Eind  <- 37.6 # DICE2023
#CDICE, Folini et al. (2023):
Ftot <- 2
Mat  <- 851
Mup  <- 628
Mlo  <- 1323
Tlo  <- 0.27
Tat  <- 1.1
H    <- 0.13

vector.ini<-list(
  ini_delc   = 0,
  ini_tildey = 0,
  ini_E      = Eind+eps_0,
  ini_Eind   = Eind,
  ini_F      = Ftot,
  ini_Mat    = Mat,                                                            
  ini_Mup    = Mup,                                                            
  ini_Mlo    = Mlo,                                                           
  ini_Tat    = Tat,                                                            
  ini_Tlo    = Tlo,
  #ini_CumE   = Eind+eps_0,
  ini_CumE   = 0,
  ini_Cumdelc= 0,
  ini_H      = H
)
remove(Eind,Ftot,Mat,Mup,Mlo,Tat,Tlo,H)


#---- Economic parameters ------------------------------------------------------
n.eta    <- 2
Phi.prep <- matrix(0,n.eta,n.eta)
n.Z <- length(vector.ini)
n.W <- n.eta + 4

param.econ<-list(
  A_bar      = NaN,                                                             #will be calibrated
  sigma_a    = NaN,                                                             #will be calibrated
  Phi        = Phi.prep,
  gamma      = 7,
  delta      = (1 - .01)^5,
  delta_K    = 0.27,                                                            
  c0         = 299                                                              #Ini consumpt. (in tn) over tstep
)
remove(Phi.prep)

#---- Climate parameters -------------------------------------------------------
#CDICE
b12     <- 0.053
b23     <- 0.0082
mateq   <- 607
mueq    <- 489
mleq    <- 1281
c1      <- 0.137 #*tstep
c3      <- 0.73
c4      <- 0.00689 #*tstep
f2co2   <- 3.45
t2co2   <- 3.25

q0  <- 135.7 #ini production DICE2023
mu0 <- 0.05

#RCP 4.5 + 6
exp.mat.2100.rcp45_6 <-1168
m0 <- exp.mat.2100.rcp45_6/mateq

param.clim<-list(
  m0         = m0, 
  a_H        = NaN,
  b_H        = NaN,
  a_N        = NaN,
  b_N        = NaN,
  kappa_N    = NaN,
  a_D        = NaN,  
  b_D        = NaN,
  mu_D       = NaN,
  mu_N       = NaN,
  eps_0      = eps_0,                                                           #DICE2023
  rho        = 0.1,                                                             #DICE2023
  mateq      = mateq,                                                 
  mueq       = mueq,                                                      
  mleq       = mleq,                                                          
  phi_0      = 0.518,                                                           #DICE2023
  phi_1      = 0.801,                                                           #DICE2023
  m_pi       = mateq,                                                           
  xi_1       = c1,                                                           
  xi_2       = c3,                                                           
  xi_3       = c4,                                                           
  nu         = t2co2,                                                       
  tau        = f2co2,                                                           
  delsigma   =-0.04,                                                            #DICE2023
  e0         = vector.ini$ini_Eind,                                             
  q0         = q0,                                                            
  mu0        = mu0,                                                           
  sigma0     = vector.ini$ini_Eind/(q0*(1-mu0)),                                    
  gsigma1    =-0.015,                                                           #DICE2023
  varphi_11  = 1-b12,                                                         
  varphi_12  = b12,                                                           
  varphi_13  = 0,                                                           
  varphi_21  = b12*mateq/mueq,                                                  
  varphi_22  = 1-b12*mateq/mueq-b23,                                          
  varphi_23  = b23,                                                          
  varphi_31  = 0,                                                           
  varphi_32  = b23*mueq/mleq,                                                 
  varphi_33  = 1-b23*mueq/mleq,                                               
  gback      = 0.05,                                                            #DICE2023
  pback      = 695,                                                             #DICE2023
  theta2     = 2.6,                                                             #DICE2023,expcost2
  b_sk       = .015/.7,# .1 or .015/.7 XXXXXXXXXXXXXXXXXXXXXXXXXXX
  tol.GN     = 10^(-6),
  eps.GN     = 10^(-5),
  mu_T       = NaN,
  mu_H       = NaN
)

param.clim$lb_dT_at <- -.1 # lower bound for dT_at

remove(b12,b23,c1,c3,c4,mateq,mueq,mleq,f2co2,
       t2co2,q0,mu0,exp.mat.2100.rcp45_6)

param<-c(param.econ,param.clim)



#---- Targets ------------------------------------------------------------------
#These are the moments the model must replicate:
target_vector <- c(
  0.95,                                                                         #consump w/ damages & T_at=2
  0.9,                                                                          #consump w/ damages & T_at=4
  0.075,                                                                        #std w/ damages & T_at=4
  .45,                                                                          #E(SLR) in 2100 & T_at=2
  .93,                                                                          #E(SLR) in 2100 & T_at=4
  (1.65-.45)/(2*qnorm(.95)),                                                    #std(SLR) in 2100 & T_at=4
  .75,                                                                          #std(T_at) in 2100
  # 36*3.667 + 446*21/1000,                                                     #E(CumN) @ T_at=2
  # 87*3.667 + 1400*21/1000,                                                    #E(CumN) @ T_at=4
  # ((141*3.667 + 2614*21/1000) - (42*3.667+836*21/1000))/2,                    #std(CumN) @ T_at=4
  # 228*3.667 + 5877*21/1000,                                                   #E(CumN) @ infinity
  56*3.667,                                                                     #E(CumN) @ T_at=2
  101*3.667,                                                                    #E(CumN) @ T_at=4
  (199*3.667 - 27*3.667)/(2*qnorm(.95)),                                        #std(CumN) @ T_at=4
  376*3.667,                                                                    #E(CumN) @ infinity
  .08,                                                                          #ini E consumption growth
  .03                                                                           #ini std dev conso growth
)

names(target_vector)<-c(
  "ECumD2",
  "ECumD4",
  "stdCumD4",
  "EH2",
  "EH4",
  "stdH4",
  "stdTat2100",
  "ECumN2",
  "ECumN4",
  "stdCumN4",
  "ECumNinf",
  "mu_c0",
  "sigma_c0"
)


#---- Names of variables in X --------------------------------------------------

names.var.X <- c("delc","y_tilde","E","E_ind","Forc","M_at",
                 "M_up","M_lo","T_at","T_lo",
                 "Cum_E","Cum_dc","H",
                 "eta_A","eta_X",
                 #"eta_E",
                 #"eta_F",
                 "D","N","dT_at","dH")


#---- log-growth rate ----------------------------------------------------------
#* must be a list of one element (mu_0) called muprice_0, and a matrix of n.X*1 (mu_1)
#* elements called muprice_1. 
#* No need to put the diagonal element (A0[13,13]=1) as a correction has been added in case.
#* Already -mu_1 in the matrix. Can think of it as mu_1.
Cum_dc1     <- matrix(0,n.Z+n.W,1)
Cum_dc1[1]  <- 1
Cum_dc      <- list(muprice_0=0,muprice_1=Cum_dc1)


#---- Define a model object ----------------------------------------------------
ini_matx<-list()
inf_matx<-list()
model<-list("parameters"=param,"vec_date"=vec_date,"tstep"=tstep,
            "MAXIT"=MAXIT,
            "n.eta"=n.eta,"n.W"=n.W,"n.Z"=n.Z,
            "Tmax"=Tmax,
            "theta0"=theta0,"horiz.2100"=horiz,"target_vector"=target_vector,
            "ini_matx"=ini_matx,"inf_matx"=inf_matx,
            "mu_c"=0, # in case we impose a growth path
            "vector.ini"=vector.ini,"Cum_dc"=Cum_dc,
            "names.var.X"=names.var.X)


#---- Perform calibration ------------------------------------------------------
print("***** starting calibration *****")
model <- solveParam4D(model)
model <- solveParam4H(model)
model <- solveParam4N(model)
model <- solveParam4c(model)
model <- solveParam4T(model)
print("***** calibration: done *****")


# gamma <- 1.45
# gamma <- 1.00001
# #gamma <- 7
# model.CRRA <- model
# model.CRRA$parameters$gamma <- gamma
# # model.CRRA$target_vector["mu_c0"] <- .04
# # model.CRRA$target_vector["sigma_c0"] <- sqrt(5*.014)
# model.CRRA <- solveParam4c(model.CRRA,indic_CRRA=TRUE)
# # model.CRRA$parameters$a_D  <- gamma * model.CRRA$parameters$a_D * 0
# # model.CRRA$parameters$b_D  <- gamma * .0018 * model.CRRA$tstep
# # model.CRRA$parameters$mu_D <- gamma * model.CRRA$parameters$mu_D * 4
# # #model$parameters$mu_T <- model$parameters$mu_T * 3
# # #model.CRRA$parameters$a_D <- 0.0001
# # #model.CRRA$parameters$b_D <- 0.0001
# # model.CRRA$parameters$b_sk <- 0
# model_CRRA_sol <- model_solve(model.CRRA,
#                               indic_mitig = T,
#                               indic_CRRA = T)
# res.scc <- scc.fct.CRRA(model_CRRA_sol,t = 0,H = 200)
# print(res.scc$SCC.CO2)
# # plot(res.scc$scc.decomp)
# # omega.ZC <- matrix(0,model_CRRA_sol$n.X,1)
# # prices.ZCRF.bonds   <- varphi(model_CRRA_sol,
# #                               omega.varphi = omega.ZC,
# #                               H = 100)
# # plot(prices.ZCRF.bonds$r.t)


tic("***** Solve Initial Model *****")
model_sol <- model_solve(model,indic_CRRA = FALSE)
#plot(model_sol$mu)
toc()

# delc <- 0
# for(i in 1:length(model_sol$omega)){
#   delc <- c(delc,model_sol$omega0[[i]][1])
# }
# plot(delc,type="l")
# compute.utility.CRRA(model_sol)


# Check plots:
H <- model$horiz.2100 + 100
#Model-implied EV
EV<-EV.fct(model_sol,H)
par(mfrow=c(2,3))
plot(model$vec_date[1:H],EV$EX$T_at,type="l")
lines(model$vec_date[1:H],EV$EX$T_atW,type="l",col="red")
lines(model$vec_date[1:H],EV$EX$T_at+2*sqrt(EV$VX$T_at),type="l",col="red",lty=2)
lines(model$vec_date[1:H],EV$EX$T_at-2*sqrt(EV$VX$T_at),type="l",col="red",lty=2)
plot(model$vec_date[1:H],EV$EX$delc,type="l")
plot(model$vec_date[1:H],EV$EX$E,type="l")
plot(model$vec_date[1:H],EV$EX$H,type="l")
lines(model$vec_date[1:H],EV$EX$H+2*sqrt(EV$VX$H),type="l",col="red",lty=2)
lines(model$vec_date[1:H],EV$EX$H-2*sqrt(EV$VX$H),type="l",col="red",lty=2)
plot(model$vec_date[1:H],EV$EX$N,type="l")

#Remove unnecessary elements:
remove(horiz,MAXIT,target_vector,vec_date,
       tstep,Tmax,n.eta,n.W,n.Z,theta0,ini_matx,
       vector.ini,param.clim,param.econ,Cum_dc,Cum_dc1,
       param,inf_matx,eps_0,H,m0,names.var.X)
