# ==============================================================================
# Figure illustrating the calibration approach
# ==============================================================================

# Plot ----
FILE = "/outputs/Figures/Figure_Damage_comparison.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=11, width=8, height=5)

par(plt=c(.2,1.0,.15,.95))

# nf <- layout(
#   matrix(c(1,1,2,3), ncol=2, byrow=TRUE), 
#   widths=c(3,1), 
#   heights=c(2,2)
# )
nf <- layout(
  matrix(c(1,2), ncol=2, byrow=TRUE), 
  widths=c(3,2), 
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
  "L. Barrage-Nordhaus",
  "L. Howard-Sterner",
  "L. Weitzman",
  "L. Barnett-Brock-Hansen",
  "L. Traeger",
  "L. Dietz-Stern"
)

all.Temp    <- seq(1.5,5,length.out=15)

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
  lines(all.Temp,All_damages[i,],col=1+i,lwd=2,type="b",pch=i)
  #points(all.Temp,All_damages[i,],col=1+i,lwd=2,pch=i)
}

# Add line for CR average damage function:
x2 <- 2
y2 <- 1 - model_sol$target_vector["ECumD2"]
x4 <- 4
y4 <- 1 - model_sol$target_vector["ECumD4"]
x0 <- 0
y0 <- y2 + (y4 - y2)/(x4 - x2)*(x0 - x2)
x10 <- 10
y10 <- y2 + (y4 - y2)/(x4 - x2)*(x10 - x2)
lines(c(x0,x10),c(y0,y10),lwd=3,col="black",lty=1)
points(c(x2,x4),c(y2,y4),pch=15,col="red",cex=1.3)

plot.new()

par(plt=c(.1,.95,.25,.85))

legend("topleft",
       legend=all.damages,
       #title="Alternative",
       lty=1,
       col = 1+(1:length(all.damages)),
       pch = 1:length(all.damages),
       #cex=1.5,
       seg.len = 3,
       lwd=2,bty = "n")

dev.off()
