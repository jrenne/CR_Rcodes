# ==============================================================================
# Figure illustrating the calibration approach
# ==============================================================================

alpha <- .032 # curvature of temperature trajectory

indic_linear_path <- FALSE # if TRUE, use linear temperature path
# (use curvature alpha otherwise)

model_sol$parameters$b_sk <- 0.01875 # 0.01875


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
                                trace_lines = FALSE,
                                lwd.CR = 3,
                                indic_add_SLR = TRUE)

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
  "L. Traeger (ACE-base)",
  "L. Traeger (ACE-HSP)",
  "L. Dietz-Stern"
)
all.colors <- c("#B22222",
                "cornsilk4",
                "#008B8B",
                "#458B00",
                "#CD5555",
                "darkorange2",
                "#838B8B",
                "darkorchid3",
                "darkolivegreen",
                "darkolivegreen",
                "darkgoldenrod3")

all.Temp    <- seq(1.5,5,length.out=15)

All_damages <- matrix(NaN,length(all.damages),length(all.Temp))

for(j in 1:length(all.Temp)){
  T2100 <- all.Temp[j]
  
  if(indic_linear_path){
    seq_of_Temp <- seq(model_sol$vector.ini$ini_Tat,T2100,length.out = h_before_2100)
  }else{
    x <- T2100
    T0 <- model_sol$vector.ini$ini_Tat
    Tinf <- (x-T0*exp(-alpha*h_before_2100))/
      (1-exp(-alpha*h_before_2100))
    seq_of_Temp <- Tinf - (Tinf - T0)*exp(-alpha * (1:h_before_2100))
  }

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
  lines(all.Temp,All_damages[i,],col=all.colors[i],lwd=2,type="b",pch=i-1)
  #points(all.Temp,All_damages[i,],col=1+i,lwd=2,pch=i)
}

points(c(2,4),
       c(1-model_sol$target_vector["ECumD2"],1-model_sol$target_vector["ECumD4"]),
       pch=15,col="red")

# Add line for CR average damage function:
# x2 <- 2
# y2 <- 1 - model_sol$target_vector["ECumD2"]
# x4 <- 4
# y4 <- 1 - model_sol$target_vector["ECumD4"]
# x0 <- 0
# y0 <- y2 + (y4 - y2)/(x4 - x2)*(x0 - x2)
# x10 <- 10
# y10 <- y2 + (y4 - y2)/(x4 - x2)*(x10 - x2)
#lines(c(x0,x10),c(y0,y10),lwd=4,col="black",lty=1)
#points(c(x2,x4),c(y2,y4),pch=15,col="red",cex=1.3)

plot.new()

par(plt=c(.1,.95,.25,.85))

legend("topleft",
       legend=c(all.damages,
                "CR (with SLR damages)",
                "CR (w/o SLR damages)",
                "Targeted mean (w/o SLR)"),
       #title="Alternative",
       lty=c(rep(NaN,length(all.damages)),1,2,NaN),
       #col = 1+(1:length(all.damages)),
       col=c(all.colors,"black","black","red"),
       pch = c(0:(length(all.damages)-1),NaN,NaN,15),
       pt.cex=1.3,
       #cex=1.5,
       seg.len = 3,
       lwd=c(rep(2,length(all.damages)),3,3),
       bty = "n")

dev.off()
