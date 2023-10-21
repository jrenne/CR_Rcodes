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


#Integers: from 0 to Tmax=100
passive.mu        <- 0  #in which period mu can reach 1?
#Binary operators: 0 = NO, 1 = YES.
indic_plots_paper <- 0  #do you run some plots of the paper? see description below
indic_tables_paper<- 0  #do you run the tables of the paper?
indic_stocks      <- 0  #do you add stocks to your analysis? see file -stocks-

#Load libraries of functions:
source("./procedure/Functions_General_v2.R")
source("./procedure/prepare.figures.R")

tic("Calibration")
source("./estimations/load_ini_model.R")
#source("./estimations/plots_check.R") 
toc()


#Updating Plots and Tables
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


plots <- 13:14
#plots <- 5

if(indic_plots_paper==1){
  source("outputs/plots_paper.R")
}

if(indic_tables_paper==1){
  print("Preparing tables (about 10 seconds)")
  source("outputs/tables_paper.R")
}

#####
#* To add a variable/parameter: see CHANGE
#* FILES:
#* load_ini_model (n.X, param, ...)
#* optimiz (Number of parameters FILTER)
#* Functions_General_v2 (model_solve, EV.fct, simul.function)
#* Functions_Optim (Param2Model, Model2Param)