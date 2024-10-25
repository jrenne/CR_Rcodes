# ==============================================================================
# Prepare tex table with solution to utility
# ==============================================================================

# Define format of figures:
nb.dec <- 3 # number of decimal numbers
Format  <- paste("%.",nb.dec,"f",sep="")
Format0 <- paste("%.",0,"f",sep="")
# Format1 <- paste("%.",1,"f",sep="")
# Format2 <- paste("%.",2,"f",sep="")
# Format3 <- paste("%.",3,"f",sep="")
# Format4 <- paste("%.",4,"f",sep="")
# Format5 <- paste("%.",5,"f",sep="")

nb.dec <- 2
Format <- paste("%3.",nb.dec,"e",sep="")

Make.entry <- function(x,Format){
  if(grepl("e-",sprintf(Format, x))){
    res <- paste("$",gsub("e-","\\\\times 10^{-",sprintf(Format, x)),"}$",sep="")
  }
  if(grepl("e+",sprintf(Format, x))){
    res <- paste("$",gsub("e+","\\\\times 10^{",sprintf(Format, x)),"}$",sep="")
  }
  plus_zero_res  <- paste("$",gsub("e+","\\\\times 10^{",sprintf(Format, 0)),"}$",sep="")
  minus_zero_res <- paste("$-",gsub("e+","\\\\times 10^{",sprintf(Format, 0)),"}$",sep="")
  if((res == minus_zero_res)|(res == plus_zero_res)){
    res <- "$-$"
  }
  return(res)
}


chosen_state_variables <- c("y_tilde",
                            "E",#"E_ind",
                            "Forc","M_at",
                            "M_up","M_lo","T_lo",
                            #"H","eta_A","eta_X","D","N",
                            "T_at"#,"HW"
                            )

matrix_names_Latex <- matrix(NaN,length(model_sol$names.var.X),2)
matrix_names_Latex[,1] <- model_sol$names.var.X
matrix_names_Latex[,2] <- c("$\\Delta c_t$","$\\tilde{y}_t$",
                            "$\\mathcal{E}_t$","$\\mathcal{E}_{Ind,t}$",
                            "$F_t$","$M_{at,t}$",
                            "$M_{up,t}$","$M_{lo,t}$",
                            "$T_{lo,t}$","$Cum_{\\mathcal{E},t}$",
                            "$Cum_{\\Delta c,t}$","$H_t$",
                            "$\\eta_{A,t}$","$\\eta_{X,t}$",
                            "$D_t$","$N_t$",
                            "$T_{at,t}$","$\\Delta H_t$")

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
                   Make.entry(model_sol$inf_matx$eta0.inf,Format),sep="")
# Short-term rate (initial period):
this.line <- paste(this.line,"&",
                   Make.entry(model_sol$eta0[1],Format),sep="")

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
                     Make.entry(model_sol$inf_matx$eta1.inf[indic],Format),sep="")
  # Short-term rate (initial period):
  this.line <- paste(this.line,"&",
                     Make.entry(model_sol$eta1[[1]][indic],Format),sep="")
  
  this.line <- paste(this.line,"\\\\",sep="")
  latex.table <- rbind(latex.table,
                       this.line)
}
  

latex.file <- paste("outputs/Tables/table_utility_solution.txt",sep="")

write(latex.table, file = latex.file)

