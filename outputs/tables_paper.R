# ==============================================================================
# LATEX TABLES
# ==============================================================================

#* Description
#* 1 Estimated parameters
#* 2 Targets + estimated param
#* 3 Initial value of State Vector
#* 4 Parameters
#* 5 Utility solution

par(plt=c(.1,.9,.1,.9))
param<-model_sol$parameters
# Prepare latex tables

make.entry <- function(x,format.nb){
  output <- paste("$",sprintf(format.nb,x),"$",sep="")
  return(output)
}

# 1 Estimated parameters
print("Preparing table of parameter estimates")
source("outputs/make_tables/make_table_Estimated_param.R")

# 2 Initial value of State Vector
print("Preparing table of initial values of state variables")
source("outputs/make_tables/make_table_IniX.R")

# 3 Parameters
print("Preparing table of parameters")
source("outputs/make_tables/make_table_param.R")

# 4 SCC sensitivity
print("Preparing table showing results of SCC sensitivity analysis (takes about 3 minutes)")
source("outputs/make_tables/make_table_SCC.R")

# 5 Solution to Utility function
print("Preparing table showing utility function specification")
source("outputs/make_tables/make_table_utility_solution.R")


