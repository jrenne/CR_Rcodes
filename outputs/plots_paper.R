#* Description
#* 0 calibration
#* 1 tibs and swaps
#* 2 digital option
#* 3 pdf temperatures
#* 4 pdf carbon concentration in atmosphere
#* 5 climate beta
#* 6 Disasters Simulations
#* 7 pdf sea level
#* 8 mu (comparison with DICE)
#* 9 Radiative forcing approximation
#*10 Constant maturity - ZCB
#*11 Cut Climate Premium
#*12 Break-even rates of inflation
#*13 Merton 1
#*14 Merton 2
#*15 Housing prices
#*16 Confidence Interval SCC and Risk premiums
#*17 Relationship between SCC and Temperature risk premium
#*18 Relationship between SCC and Temperature risk premium, Lemoine approach
#*19 Illustrations of gamma-zero distribution
#*20 Comparison of damages
#*21 IRF 1 Gt Carbon
#*22 Comparison of SCC with ACE model
#*23 Comparison of risk-free yield curves with alternative approaches

# Maturity:
H <-model_sol$horiz.2100                                                        #maturity 2100
EV<-EV.fct(model_sol,H)
#Pricing
omega_ZCB  <- model_sol$omega_ZCB
omega_T.at <- model_sol$omega_T.at

a <- omega_T.at                                                                 #Options related to T_at in X
b <- 2                                                                          #Options pay when a'X < b

# Miscellaneous
n.date<-model_sol$Tmax

# Quantiles considered for pdf charts:
vector.of.CI <- c(0.5,0.8,0.9,0.95)

# Colors:
P.col.line <- brocolors("crayons")["Teal Blue"]
P.col<-adjustcolor( P.col.line, alpha.f = 0.15)
Q.col.line <- brocolors("crayons")["Mango Tango"]
Q.col<-adjustcolor( Q.col.line, alpha.f = 0.15)

# Nb of values used for refined grid when computing Conf Intervals:
nb.values.variable <- 400



#********************************0*********************************************#
#CALIBRATION
if(is.element(0,plots)){
  print("Preparing Calibration plot")
  source("outputs/make_figures/make_figure_calibration.R")
}

#********************************1*********************************************#
#TIBS/SWAPS
if(is.element(1,plots)){
  print("Preparing TS plot")
  source("outputs/make_figures/make_figure_TS.R")
}

#********************************2*********************************************#
#Digital Option
if(is.element(2,plots)){
  print("Preparing Option plot")
  source("outputs/make_figures/make_figure_options.R")
}

#********************************3*********************************************#
#Temperatures pdf + RCP
if(is.element(3,plots)){
  print("Preparing Temperature distri plot")
  source("outputs/make_figures/make_figure_Tpdf.R")
}

#********************************4*********************************************#
#Carbon Concentration Atmosphere pdf + RCP
if(is.element(4,plots)){
  print("Preparing Atm carbon concentration distri plot")
  source("outputs/make_figures/make_figure_Mpdf.R")
}

#********************************5*********************************************#
#Sensitivity to mu_D
if(is.element(5,plots)){
  print("Preparing sensitivity-to-mu_D plot (takes about 30 sec)")
  source("outputs/make_figures/make_figure_sensitiv_muD.R")
  source("outputs/make_figures/make_figure_SCC.R")
}

#********************************6*********************************************#
#Disasters Simulations
if(is.element(6,plots)){
  print("Preparing D simulation plot")
  source("outputs/make_figures/make_figure_Dsimul.R")
}

#********************************7*********************************************#
#Global Sea Level pdf + RCP
if(is.element(7,plots)){
  print("Preparing SLR distri plot")
  source("outputs/make_figures/make_figure_Hpdf.R")
}

#********************************8*********************************************#
#Mitigation vs DICE2016
if(is.element(8,plots)){
  print("Preparing figure mu plot (takes about 30 sec)")
  source("outputs/make_figures/make_figure_mu.R")
}

#*******************************9*********************************************#
#Radiative Forcings Approximation
if(is.element(9,plots)){
  print("Preparing radiative forcing approx plot")
  source("outputs/make_figures/make_figure_RFapprox.R")
}

#*******************************10*********************************************#
#Constant maturity for ZCB
if(is.element(10,plots)){
  print("Preparing ZCB plot")
  source("outputs/make_figures/make_figure_ConstantMaturityZCB.R")
}


#*******************************11*********************************************#
#Cut in Climate Premium
if(is.element(11,plots)){
  print("Preparing plot of sensitivitiy of Temp risk premium to mu_D")
  source("outputs/make_figures/make_figure_cut_CP_muD.R")
}

#*******************************12*********************************************#
#Break-even rates of inflation
if(is.element(12,plots)){
  print("Preparing figure BEIR term structures")
  source("outputs/make_figures/make_figure_breakeveninflation.R")
}

#*******************************13-14******************************************#
#Merton model
if(is.element(13,plots)){
  print("Preparing Merton-model figure (1/2)")
  source("outputs/make_figures/make_figure_Merton.R")
}
if(is.element(14,plots)){
  print("Preparing Merton-model figure (2/2)")
  source("outputs/make_figures/make_figure_Merton2.R")
}

#*******************************15*********************************************#
#Housing prices
if(is.element(15,plots)){
  print("Preparing figure Housing")
  source("outputs/make_figures/make_figure_Housing.R")
}

#*******************************16*********************************************#
#Confidence intervals for SCC and risk premiums
if(is.element(16,plots)){
  print("Preparing figure Confidence int. for SCC and Risk premiums (takes a few minutes)")
  source("outputs/make_figures/make_figure_ConfInt_RP.R",
         encoding = 'ISO8859-1')
}

#*******************************17*********************************************#
#Relationship between SCC and Temperature risk premium
if(is.element(17,plots)){
  print("Preparing figure Relationship SCC and Temp. risk premium (takes a few minutes)")
  source("outputs/make_figures/make_figure_SCC_vs_TempRP.R",
         encoding = 'ISO8859-1')
}

#*******************************18*********************************************#
#Relationship between SCC and Temperature risk premium, Lemoine's approach
if(is.element(18,plots)){
  print("Preparing figure Relationship SCC and Temp. risk premium, Lemoine's approach")
  source("outputs/make_figures/make_figure_SCC_vs_TempRP_Lemoine.R",
         encoding = 'ISO8859-1')
}

#*******************************19*********************************************#
#illustrations of the Gamma-zero distribution
if(is.element(19,plots)){
  print("Preparing figure illustrating Gamma-zero distribution")
  source("outputs/make_figures/make_figure_gamma0_distri.R",
         encoding = 'ISO8859-1')
}

#*******************************20*********************************************#
#comparison of damage specifications
if(is.element(20,plots)){
  print("Preparing figure illustrating Gamma-zero distribution")
  source("outputs/make_figures/make_figure_Damage_comparison.R",
         encoding = 'ISO8859-1')
}

#*******************************21*********************************************#
#IRF 1 Gt Carbon
if(is.element(21,plots)){
  print("Preparing figure showing temperature dynamic effect of 1Gt carbon pulse")
  source("outputs/make_figures/make_figure_IRF100Gt.R",
         encoding = 'ISO8859-1')
}

#*******************************22*********************************************#
#Comparison of SCC with ACE model
if(is.element(22,plots)){
  print("Preparing figure comparing SCC with ACE model (takes about 3 minutes)")
  source("outputs/make_figures/make_figure_compar_SCC_Traeger.R",
         encoding = 'ISO8859-1')
}

#*******************************23*********************************************#
#Comparison of risk-free yield curves with alternative approaches
if(is.element(23,plots)){
  print("Preparing figure comparing risk-free yield curves (takes about 20 seconds)")
  source("outputs/make_figures/make_figure_YC_RF.R",
         encoding = 'ISO8859-1')
}


