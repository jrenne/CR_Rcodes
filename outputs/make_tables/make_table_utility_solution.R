# ==============================================================================
# TABLE 2. Utility solution
# table_utility_solution.txt
# ==============================================================================

# Define format of figures:
nb.dec <- 3 # number of decimal numbers
Format  <- paste("%.",nb.dec,"f",sep="")
Format0 <- paste("%.",0,"f",sep="")

nb.dec <- 3
Format <- paste("%3.",nb.dec,"e",sep="")

Make.entry <- function(x,Format){
  if(grepl("+00",sprintf(Format, x))){
    res <- paste("$",sprintf(paste("%0.",nb.dec,"f",sep=""),x),"$",sep="")
  }else{
    if(grepl("e-",sprintf(Format, x))){
      res <- paste("$",gsub("e-","\\\\times 10^{-",sprintf(Format, x)),"}$",sep="")
    }
    if(grepl("e+",sprintf(Format, x))){
      res <- paste("$",gsub("e+","\\\\times 10^{",sprintf(Format, x)),"}$",sep="")
    }
  }
  if((res == paste("$",sprintf(paste("%0.",nb.dec,"f",sep=""),0),"$",sep=""))|
     res == paste("$-",sprintf(paste("%0.",nb.dec,"f",sep=""),0),"$",sep="")){
    res <- "$-$"
  }
  return(res)
}


chosen_state_variables <- c("y_tilde",
                            "E",
                            "Forc","M_at",
                            "M_up","M_lo","T_at","T_lo"
                            )

matrix_names_Latex <- matrix(NaN,length(model_sol$names.var.X),2)
matrix_names_Latex[,1] <- model_sol$names.var.X
matrix_names_Latex[,2] <- c("$\\Delta c_t$","$\\tilde{y}_t$",
                            "$\\mathcal{E}_t$",
                            "$F_t$","$M_{AT,t}$",
                            "$M_{UP,t}$","$M_{LO,t}$",
                            "$T_{AT,t}$","$T_{LO,t}$",
                            "$Cum_{\\mathcal{E},t}$",
                            "$Cum_{\\Delta c,t}$","$H_t$",
                            "$\\eta_{A,t}$","$\\eta_{X,t}$",
                            "$D_t$","$N_t$",
                            "$T_{AT,t}$","$\\Delta H_t$")

# Compute utility solution for initial period:
mu_u.1 <- mu_u.t.fct(model_sol,model_sol$Tmax)

# Initialize table:
latex.table <- NULL

# Constants: -------------------------------------------------------------------
this.line <- "$\\mu_{u,0,t}$/$\\mu_{r,0,t}$"
# Utility (after Tmax):
this.line <- paste(this.line,"&",
                   Make.entry(model_sol$mu_u0,Format),sep="")
# Utility (initial period):
this.line <- paste(this.line,"&",
                   Make.entry(mu_u.1$mu_u0.t,Format),sep="")
# Short-term rate (after Tmax):
this.line <- paste(this.line,"&",
                   Make.entry(model_sol$inf_matx$eta0.inf/model_sol$tstep,Format),sep="")
# Short-term rate (initial period):
this.line <- paste(this.line,"&",
                   Make.entry(model_sol$eta0[1]/model_sol$tstep,Format),sep="")

this.line <- paste(this.line,"\\\\",sep="")
latex.table <- rbind(latex.table,
                     this.line,
                     "\\hline")

for(i in 1:length(chosen_state_variables)){
  
  indic <- which(chosen_state_variables[i]==matrix_names_Latex[,1])
  this.line <- matrix_names_Latex[indic,2]
  
  # Utility (after Tmax):
  this.line <- paste(this.line,"&",
                     Make.entry(model_sol$mu_u1[indic],Format),sep="")
  # Utility (initial period):
  this.line <- paste(this.line,"&",
                     Make.entry(mu_u.1$mu_u1.t[indic],Format),sep="")
  # Short-term rate (after Tmax):
  this.line <- paste(this.line,"&",
                     Make.entry(model_sol$inf_matx$eta1.inf[indic]/model_sol$tstep,Format),sep="")
  # Short-term rate (initial period):
  this.line <- paste(this.line,"&",
                     Make.entry(model_sol$eta1[[1]][indic]/model_sol$tstep,Format),sep="")
  
  this.line <- paste(this.line,"\\\\",sep="")
  latex.table <- rbind(latex.table,
                       this.line)
}
  

latex.file <- paste("outputs/Tables/table_utility_solution.txt",sep="")

write(latex.table, file = latex.file)

