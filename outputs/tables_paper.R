# ==============================================================================
# LATEX TABLES
# ==============================================================================

#* Description
#* 1 Estimated parameters
#* 2 Targets + estimated param
#* 3 Initial value of State Vector
#* 4 Parameters


par(plt=c(.1,.9,.1,.9))
param<-model_sol$parameters
# Prepare latex tables

make.entry <- function(x,format.nb){
  output <- paste("$",sprintf(format.nb,x),"$",sep="")
  return(output)
}

# 1 Estimated parameters
source("outputs/make_tables/make_table_Estimated_param.R")

# 2 Initial value of State Vector
source("outputs/make_tables/make_table_IniX.R")

# 3 Parameters
source("outputs/make_tables/make_table_param.R")

# 4 SCC sensitivity
source("outputs/make_tables/make_table_SCC.R")


