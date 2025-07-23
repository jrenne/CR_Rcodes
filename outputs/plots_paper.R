#* Figures description:
#*  2-- FIGURE 2. Atmospheric temperature response to a carbon pulse 
#*      make_figure_IRF1Gt.R
#*  3-- FIGURE 3. Damage function 
#*      make_figure_Damage_comparison.R
#*  4-- FIGURE 4. Conditional distribution of future temperatures 
#*      make_figure_Tpdf.R
#*  5-- FIGURE 5. Conditional distribution of future global sea level 
#*      make_figure_Hpdf.R
#*  6-- FIGURE 6. The term structure of real rates
#*      make_figure_YC_RF.R
#*  7-- FIGURE 7. The term structures of interest rates 
#*      make_figure_breakeveninflation.R
#*  8-- FIGURE 8. Price of digital options, with contributions of risk premiums 
#*      make_figure_options.R
#*  9-- FIGURE 9. Effect of \mu on the conditional distribution x_t|y_t \sim \gamma0 (y_t/\mu, \mu)
#*      make_figure_gamma0_distri.R 
#* 10-- FIGURE 10. Distribution of cumulated damages 
#*      make_figure_gamma0_distri.R
#* 11-- FIGURE 11. From carbon concentrations to atmospheric temperature 
#*      make_figure_RCP_to_TAT.R
#* 101-- FIGURE S.1. Model calibration 
#*      make_figure_calibration.R
#* 102-- FIGURE S.2. Mitigation rate 
#*      make_figure_mu.R
#* 103-- FIGURE S.3. Social Cost of Carbon and temperature risk premiums
#*      make_figure_SCC_vs_TempRP.R
#* 104-- FIGURE S.4. Comparison of temperature risk premiums
#*      make_figure_YC_RF.R
#* 105-- FIGURE S.5. Merton model 
#*      make_figure_Merton.R
#* -----------------------------------------------------------------------------

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

#--  6 + 103
# FIGURE 6. The term structure of real rates 
# + FIGURE S.4. Comparison of temperature risk premiums
if(is.element(6|104,plots)){
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

#-- 101
# FIGURE S.1. Model calibration
if(is.element(101,plots)){
  print("Preparing Calibration plot")
  source("outputs/make_figures/make_figure_calibration.R")
}

#-- 102
# FIGURE S.2. Mitigation rate
if(is.element(102,plots)){
  print("Preparing figure mu plot (takes about 30 sec)")
  source("outputs/make_figures/make_figure_mu.R")
}

#-- 103
# FIGURE S.4. Social Cost of Carbon and temperature risk premiums
if(is.element(103,plots)){
  print("Preparing figure Relationship SCC and Temp. risk premium (takes a few minutes)")
  source("outputs/make_figures/make_figure_SCC_vs_TempRP.R",
         encoding = 'ISO8859-1')
}

#-- 105
# FIGURE S.5. Merton model
if(is.element(105,plots)){
  print("Preparing Merton-model figure")
  source("outputs/make_figures/make_figure_Merton.R")
}
