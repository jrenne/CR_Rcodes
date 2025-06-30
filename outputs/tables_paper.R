# ==============================================================================
# LATEX TABLES
# ==============================================================================

#* Description:
#* TABLE 1. Model-specific parameters 
#* TABLE 2. Utility solution
#* TABLE 5. Social Cost of Carbon, in US $ per ton of CO2
#* TABLE 7. Calibrated parameters
#* TABLE VI.1. 2100 Temperature risk premium
#* TABLE VI.2. Long-term interest rate (maturity: 2100)

par(plt=c(.1,.9,.1,.9))
param<-model_sol$parameters
# Prepare latex tables

make.entry <- function(x,format.nb){
  output <- paste("$",sprintf(format.nb,x),"$",sep="")
  return(output)
}

# TABLE 1. Model-specific parameters
print("Preparing table of parameter estimates")
source("outputs/make_tables/make_table_Estimated_param.R")

# TABLE 2. Utility solution
print("Preparing table showing utility function specification")
source("outputs/make_tables/make_table_utility_solution.R")

# TABLE 5. Social Cost of Carbon, in US $ per ton of CO2
# TABLE VI.1. 2100 Temperature risk premium
# TABLE VI.2. Long-term interest rate (maturity: 2100)
print("Preparing table showing results of SCC sensitivity analysis (takes about 3 minutes)")
source("outputs/make_tables/make_table_SCC.R")

# TABLE 7. Calibrated parameters
print("Preparing table of parameters")
source("outputs/make_tables/make_table_param.R")








