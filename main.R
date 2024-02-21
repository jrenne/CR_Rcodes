# ==============================================================================
# An Analytical Framework to Price Long-Dated Climate-Exposed Assets
# ------------------------------------------------------------------------------
# Pauline Chikhani and Jean-Paul Renne
# This version: February 2024
# ==============================================================================


#clear environment
rm(list=ls(all=T)) 
library(tictoc)
library(parallel)
library(doParallel)
library(mgcv)
library(colorspace)
library(broman)
library(optimx)
library(MASS)
library(expm)


#Binary operators: 0 = NO, 1 = YES.
indic_plots_paper  <- 1  #do you run some plots of the paper? see description below
indic_tables_paper <- 1  #do you run the tables of the paper?

# For scripts using parallel computing:
number.of.cores <- 8

#Load libraries of functions:
source("procedures/functions_general.R")
source("procedures/functions_figures.R")

tic("Calibration")
source("estimations/load_ini_model.R")
#source("estimations/plots_check.R") 
toc()


#Updating Plots and Tables -----------------------------------------------------

#* Description:
#* 0 calibration
#* 1 TIBS and swaps
#* 2 digital option
#* 3 pdf temperatures
#* 4 pdf carbon concentration in atmosphere
#* 5 climate beta
#* 6 Disasters Simulations
#* 7 pdf sea level
#* 8 mu (comparison with DICE)
#* 9 Radiative forcing approximation
#*10 Constant maturity - ZCB
#*11 Sensitivity of climate Premium
#*12 Break-even rates of inflation
#*13 Merton 1
#*14 Merton 2
#*15 Housing prices
#*16 Confidence Interval SCC and Risk premiums
#*17 Relationship between SCC and Temperature risk premium
#*18 Relationship between SCC and Temperature risk premium, alternative models

plots <- 0:18

if(indic_plots_paper==1){
  source("outputs/plots_paper.R")
}

if(indic_tables_paper==1){
  print("Preparing tables (couple of minutes)")
  source("outputs/tables_paper.R")
}

