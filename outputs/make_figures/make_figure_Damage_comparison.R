# ==============================================================================
# Figure illustrating the calibration approach
# ==============================================================================

# Plot ----
FILE = "/outputs/Figures/Figure_Damage_comparison.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=11, width=9, height=6)

par(plt=c(.2,1.0,.15,.95))

# nf <- layout(
#   matrix(c(1,1,2,3), ncol=2, byrow=TRUE), 
#   widths=c(3,1), 
#   heights=c(2,2)
# )
nf <- layout(
  matrix(c(1,2), ncol=2, byrow=TRUE), 
  widths=c(2,1), 
  heights=c(1)
)

# Damages ----------------------------------------------------------------------
make_figure_calibration_Damages(model_sol,
                                main.title = "",
                                trace_lines = FALSE)


h_before_2100 <- which(model_sol$vec_date==2100) - 1

all.damages <- c(
  "G. Lemoine",
  "G. Dell-Jones-Olken",
  "G. Bansal-Kiku-Ochoa",
  "L. Nordhaus-Sztorc",
  "L. Weitzman",
  "L. Barnett-Brock-Hansen",
  "L. Traeger",
  "L. Dietz-Stern"
)

all.Temp    <- seq(1.5,5,length.out=30)

All_damages <- matrix(NaN,length(all.damages),length(all.Temp))

for(j in 1:length(all.Temp)){
  T2100 <- all.Temp[j]
  seq_of_Temp <- seq(model_sol$vector.ini$ini_Tat,T2100,length.out = h_before_2100)
  
  for(i in 1:length(all.damages)){
    
    damage.type <- all.damages[i]
    
    if(stringr::str_sub(damage.type,1,1)=="G"){
      cumD <- prod(1 + compute_alternative_damage(T = seq_of_Temp,
                                                  damage.type,
                                                  time_step = model_sol$tstep))-1  
    }else{
      cumD <- 1 - 
        compute_alternative_damage(T = seq_of_Temp,damage.type)[h_before_2100]
    }
    
    All_damages[i,j] <- cumD
  }
}

for(i in 1:length(all.damages)){
  lines(all.Temp,All_damages[i,],col=1+i,lwd=2)
}

plot.new()

par(plt=c(.1,.95,.25,.85))

legend("topleft",
       legend=all.damages,
       #title="Alternative",
       lty=1,
       col = 1+(1:length(all.damages)),
       #cex=1.5,
       lwd=2,bty = "n")

dev.off()
