#* Figures description:
#* -----------------------------------------------------------------------------
#*  2-- FIGURE 2. Atmospheric temperature response to a carbon pulse
#*  3-- FIGURE 3. Damage function
#*  4-- FIGURE 4. Conditional distribution of future temperatures 
#*  5-- FIGURE 5. Conditional distribution of future global sea level 
#*  6-- FIGURE 6. The term structure of real rates
#*  7-- FIGURE 7. The term structures of interest rates 
#*  8-- FIGURE 8. Price of digital options, with contributions of risk premiums 
#*  9-- FIGURE 9. Effect of \mu on the conditional distribution x_t|y_t \sim \gamma0 (y_t/\mu, \mu) 
#* 10-- FIGURE 10. Distribution of cumulated damages 
#* 11-- FIGURE 11. From carbon concentrations to atmospheric temperature 
#* 31-- FIGURE III.1. Model calibration 
#* 51-- FIGURE V.1. Mitigation rate 
#* 52-- FIGURE V.2. Social Cost of Carbon and temperature risk premiums
#* 53-- FIGURE V.3. Comparison of temperature risk premiums 
#* 54-- FIGURE V.4. Merton model 
#* -----------------------------------------------------------------------------
 
 
#* Miscellaneous

#* Figure 1*: From carbon concentrations to atmospheric temperature
# print("Preparing figure showing effect of linearization on 2300 temp. distri")
# source("outputs/make_figures/make_figure_distrTAT_4_simul_Mat.R",
#        encoding = 'ISO8859-1')

#* Figure 2*: Simulated damages
# source("outputs/make_figures/make_figure_simul_cumDamages.R",
#        encoding = 'ISO8859-1')

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

#--  2
# FIGURE 2. Atmospheric temperature response to a carbon pulse
if(is.element(2,plots)){
  print("Preparing figure showing temperature dynamic effect of 1Gt carbon pulse")
  source("outputs/make_figures/make_figure_IRF1Gt.R",
         encoding = 'ISO8859-1')
}

#--  3
# FIGURE 3. Damage function
if(is.element(3,plots)){
  print("Preparing figure illustrating Gamma-zero distribution")
  source("outputs/make_figures/make_figure_Damage_comparison.R",
         encoding = 'ISO8859-1')
}

#--  4
# FIGURE 4. Conditional distribution of future temperatures
if(is.element(4,plots)){
  print("Preparing Temperature distri plot")
  source("outputs/make_figures/make_figure_Tpdf.R")
}

#--  5
# FIGURE 5. Conditional distribution of future global sea level
if(is.element(5,plots)){
  print("Preparing SLR distri plot")
  source("outputs/make_figures/make_figure_Hpdf.R")
}

#--  6 + 53
# FIGURE 6. The term structure of real rates 
# + FIGURE V.3. Comparison of temperature risk premiums
if(is.element(6|53,plots)){
  print("Preparing figure comparing risk-free yield curves (takes about 20 seconds)")
  source("outputs/make_figures/make_figure_YC_RF.R",
         encoding = 'ISO8859-1')
}

#--  7
# FIGURE 7. The term structures of interest rates
if(is.element(7,plots)){
  print("Preparing figure BEIR term structures")
  source("outputs/make_figures/make_figure_breakeveninflation.R")
}

#--  8
# FIGURE 8. Price of digital options, with contributions of risk premiums
if(is.element(8,plots)){
  print("Preparing Option plot")
  source("outputs/make_figures/make_figure_options.R")
}

#--  9 + 10
# FIGURE 9. Effect of \mu on the conditional distribution x_t|y_t \sim \gamma0 (y_t/\mu, \mu)
# + FIGURE 10. Distribution of cumulated damages
if(is.element(9|10,plots)){
  print("Preparing figure illustrating Gamma-zero distribution")
  source("outputs/make_figures/make_figure_gamma0_distri.R",
         encoding = 'ISO8859-1')
}

#--  11
# FIGURE 11. From carbon concentrations to atmospheric temperature
if(is.element(11,plots)){
  print("Preparing figure showing response of temperature conditional on MAT paths")
  source("outputs/make_figures/make_figure_RCP_to_TAT.R",
         encoding = 'ISO8859-1')
}

#-- 31
# FIGURE III.1. Model calibration
if(is.element(31,plots)){
  print("Preparing Calibration plot")
  source("outputs/make_figures/make_figure_calibration.R")
}

#-- 51
# FIGURE V.1. Mitigation rate
if(is.element(51,plots)){
  print("Preparing figure mu plot (takes about 30 sec)")
  source("outputs/make_figures/make_figure_mu.R")
}

#-- 52
# FIGURE V.2. Social Cost of Carbon and temperature risk premiums
if(is.element(52,plots)){
  print("Preparing figure Relationship SCC and Temp. risk premium (takes a few minutes)")
  source("outputs/make_figures/make_figure_SCC_vs_TempRP.R",
         encoding = 'ISO8859-1')
}

#-- 54
# FIGURE V.4. Merton model
if(is.element(54,plots)){
  print("Preparing Merton-model figure")
  source("outputs/make_figures/make_figure_Merton.R")
}
