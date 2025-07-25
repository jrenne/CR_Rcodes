# ==============================================================================
# An Analytical Framework to Price Long-Dated Climate-Exposed Assets
# ------------------------------------------------------------------------------
# Pauline Chikhani and Jean-Paul Renne (jean-paul.renne@unil.ch)
# This version: July 2025
# ==============================================================================

# Clear environment ------------------------------------------------------------
rm(list=ls(all=T))

# Set working directory --------------------------------------------------------
setwd("~/Dropbox/Research/TIBs/CR_Rcodes")

# Load libraries ---------------------------------------------------------------
library(parallel)
library(doParallel)
library(colorspace)
library(broman)
library(optimx)
library(MASS)
library(expm)

# Determine whether tables and figures are generated ---------------------------
#Binary operators: 0 = NO, 1 = YES.
indic_plots_paper  <- 1 #produce paper's plots? see description below
indic_tables_paper <- 1 #produce paper's table? see description below

# Set number of cores used for parallel computing (sensitivity analysis) -------
number.of.cores <- 8

# Load libraries of functions --------------------------------------------------
source("procedures/functions_general.R")
source("procedures/functions_figures.R")
source("procedures/functions_other_models.R")

# Calibrate baseline model -----------------------------------------------------
source("estimations/load_ini_model.R")


# Generate Figures and Tables --------------------------------------------------

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

plots <- 1:105 # this selects all possible figures

if(indic_plots_paper==1){
  source("outputs/plots_paper.R")
}

#* Tables description:
#* TABLE 1. Model-specific parameters 
#*          make_table_Estimated_Param.R
#* TABLE 2. Utility solution
#*          make_table_utility_solution.R
#* TABLE 5. Social Cost of Carbon, in US $ per ton of CO2   
#*          make_table_SCC.R
#* TABLE 7. Calibrated parameters
#*          make_table_param.R
#*          
#* TABLE S.1. 2100 Temperature risk premium
#*            make_table_SCC.R
#* TABLE S.2. Long-term interest rate (maturity: 2100)
#*            make_table_SCC.R

if(indic_tables_paper==1){
  print("Preparing tables (couple of minutes)")
  source("outputs/tables_paper.R")
}

