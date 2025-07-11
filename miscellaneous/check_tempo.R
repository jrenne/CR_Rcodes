
# ==============================================================================
# FIGURE for reply (second round)
# Figure_RCP_to_TAT.pdf
# ==============================================================================


indic_print_values <- 0

indic_Mat <- which(model_sol$names.var.X=="M_at")

RCP_MAGICC <- read.csv("data/RCP_Mat_MAGICC.csv", header=FALSE)
RCP_ACE    <- read.csv("data/RCP_Mat_ACE.csv", header=FALSE)

H2500 <- model$horiz.2100 + (2500 - 2100)/model_sol$tstep
years <- seq(model_sol$vec_date[1],2300,by=model_sol$tstep)

indicators <- which(RCP_MAGICC$V1 %in% years) # to find appropriate years

MAT_RCP30 <- RCP_MAGICC$V2[indicators]
MAT_RCP45 <- RCP_MAGICC$V3[indicators]
MAT_RCP60 <- RCP_MAGICC$V4[indicators]
MAT_RCP85 <- RCP_MAGICC$V5[indicators]

T_RCP30 <- RCP_MAGICC$V6[indicators]
T_RCP45 <- RCP_MAGICC$V7[indicators]
T_RCP60 <- RCP_MAGICC$V8[indicators]
T_RCP85 <- RCP_MAGICC$V9[indicators]

RF_RCP30 <- RCP_MAGICC$V14[indicators]
RF_RCP45 <- RCP_MAGICC$V15[indicators]
RF_RCP60 <- RCP_MAGICC$V16[indicators]
RF_RCP85 <- RCP_MAGICC$V17[indicators]

T_ACE_RCP30 <- RCP_ACE$V2
T_ACE_RCP45 <- RCP_ACE$V3
T_ACE_RCP60 <- RCP_ACE$V4
T_ACE_RCP85 <- RCP_ACE$V5

# Compute T_AT with ACE under CR -----------------------------------------------
# Compute E(M_AT) under CR baseline model:
EV <- EV.fct(model_sol)
MAT_CR <- c(model_sol$vector.ini$ini_Mat,
            EV$EX$M_at)
years_CR <- c(2020,EV$date)
indicators_CR <- which(years_CR %in% years) # to find appropriate years
MAT_CR <- MAT_CR[indicators_CR]
# Compute weight in RCP45 (versus RCP60) to approximate MAT in CR with
# linear combination of RCP45 and RCP60:

Temp.ini <- matrix(c(model_sol$vector.ini$ini_Tat,
                     model_sol$vector.ini$ini_Tlo,0),3,2)

# Determine Gt process in Traeger (2023), using his eq.(6):
eta  <- 3.8
Mpre <- 596.4 # from Matlab codes
RF <- cbind(RF_RCP45,RF_RCP60,RF_RCP85)
M1 <- cbind(MAT_RCP45,MAT_RCP60,MAT_RCP85)
G_ACE <- exp(RF*log(2)/eta)*Mpre - M1
G_ACE_RCP45 <- G_ACE[,1] 
G_ACE_RCP60 <- G_ACE[,2]
G_ACE_RCP85 <- G_ACE[,3]

# Determine RF using Traeger (2023) eq.(6):
RF_CR_4ACE    <- eta/log(2)*log((MAT_CR    + G_ACE_RCP60)/Mpre)
RF_RCP45_4ACE <- eta/log(2)*log((MAT_RCP45 + G_ACE_RCP45)/Mpre)
RF_RCP60_4ACE <- eta/log(2)*log((MAT_RCP60 + G_ACE_RCP60)/Mpre)
RF_RCP85_4ACE <- eta/log(2)*log((MAT_RCP85 + G_ACE_RCP85)/Mpre)

forcing_CR    <- matrix(1,2,1) %*% t(exp(log(2) / eta * RF_CR_4ACE))
forcing_RCP45 <- matrix(1,2,1) %*% t(exp(log(2) / eta * RF_RCP45_4ACE))
forcing_RCP60 <- matrix(1,2,1) %*% t(exp(log(2) / eta * RF_RCP60_4ACE))
forcing_RCP85 <- matrix(1,2,1) %*% t(exp(log(2) / eta * RF_RCP85_4ACE))

res <- TempSimulation_ACE(Temp.ini, forcing_CR)
T_ACE_CR <- res[1,1,]
res <- TempSimulation_ACE(Temp.ini, forcing_RCP45)
T_ACE_RCP45 <- res[1,1,]
res <- TempSimulation_ACE(Temp.ini, forcing_RCP60)
T_ACE_RCP60 <- res[1,1,]
res <- TempSimulation_ACE(Temp.ini, forcing_RCP85)
T_ACE_RCP85 <- res[1,1,]


# ------------------------------------------------------------------------------
# Plot ----
FILE = "/outputs/Figures/Figure_RCP_to_TAT_TEMPO.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=12, width=10, height=10)

par(plt=c(.06,.99,.16,.86))

par(mfrow=c(2,1))

for(case in c(1,2)){
  
  
  # Grid of values of M_AT (for first and second panel):
  if(case == 1){
    gridMat<-seq(900,3200,by=1)
    
    main.t <- expression(paste("AS SHOWN IN THE PAPER"))
      
  }else{
    gridMat<-seq(900,1700,by=1)
    
    main.t <- expression(paste("FOCUSSING ON LOWER CONCENTRATIONS"))
    
  }
  
  nlinF <- matrix(0,length(gridMat),1)
  nlinF[1:length(nlinF)] <- model_sol$parameters$tau/log(2)*
    log(gridMat[1:length(gridMat)]/model_sol$parameters$m_pi)
  
  vector.of.dates <- c(2070,2100,2150)
  colors <- c("#A2CD5A",
              "#CD3333",
              "#7AC5CD")
  
  legend_names <- NULL
  
  specific_MAT_RCP45 <- MAT_RCP45[which(years %in% vector.of.dates)]
  specific_MAT_RCP60 <- MAT_RCP60[which(years %in% vector.of.dates)]
  specific_MAT_RCP85 <- MAT_RCP85[which(years %in% vector.of.dates)]
  
  matrix_MAT_RCP <- rbind(specific_MAT_RCP45,
                          specific_MAT_RCP60,
                          specific_MAT_RCP85)
  
  ylim <- c(1,10)
  
  plot(0, 0,col="white",ylim=ylim,xlim=range(gridMat),type="l",
       las=1,xlab=expression(paste("M"[AT]," (GtC)")),
       ylab="",
       main=main.t)
  
  title(ylab=expression(paste("FCO"[2]," (in Wm-2)")), line=1.5,
        cex.lab=1)
  
  
  # Distributions of M_AT for CR model -------------------------------------------
  # For Fourier-transform computation:
  x <- exp(seq(-10,10,length.out = 5000)) # grid for Fourier
  indic.Mat <- which(model_sol$names.var.X=="M_at")
  for(iii in 1:length(vector.of.dates)){
    H_considered <- which(model_sol$vec_date==vector.of.dates[iii])-1
    
    Mat.distr.considered <- fourier(model_sol,x,gridMat,H_considered,indic.Mat)
    Mat.pdf.considered   <- diff(Mat.distr.considered)
    
    if(iii==1){
      # determine multiplication factor (for chart)
      max.pdf <- max(Mat.pdf.considered)
      multi.factor <- .8 * (ylim[2] - ylim[1])/max.pdf
    }
    
    Mat.pdf.considered <- ylim[1] + multi.factor * Mat.pdf.considered
    
    gridF    <- seq(0.00,8,by=.002)
    gridMat_1<- gridMat[-1]
    nF       <- length(gridF)
    nMat     <- length(gridMat_1)
    mx.F     <- t(matrix(gridF,nF,nMat))
    mx.Mat   <-   matrix(gridMat_1,nMat,nF)
    
    colors_tranp <- paste(colors,"77",sep="")
    
    polygon(c(gridMat[-1],rev(gridMat[-1])),
            c(Mat.pdf.considered,rep(ylim[1],length(Mat.pdf.considered))),
            col=colors_tranp[iii],border = NA)
    
    indic.considered.in.yearsvector <- which(years==vector.of.dates[iii])
    indic_quantile <- which.min((MAT_RCP85[indic.considered.in.yearsvector] - gridMat)^2)
    quantile_RCP85_considered <- Mat.distr.considered[indic_quantile]
    
    indic_quantile <- which.min((MAT_RCP60[indic.considered.in.yearsvector]+500 - gridMat)^2)
    quantile_RCP60plus500_considered <- Mat.distr.considered[indic_quantile]
    
  }
  # ------------------------------------------------------------------------------
  
  
  
  lines(gridMat,nlinF,
        col="black",lwd=2,lty=2)
  
  for(iii in 1:length(vector.of.dates)){
    
    legend_names <- c(legend_names,paste("Linearized, for year ",vector.of.dates[iii],sep=""))
    
    horiz <- vector.of.dates[iii]
    y <- which(horiz == model_sol$vec_date)
    m0 <- model_sol$parameters$m0[y]
    
    gridMat_reduced <- seq(gridMat[1],specific_MAT_RCP85[iii],by=10)
    
    linF  <-matrix(0,length(gridMat_reduced),1)
    linF[1:length(linF)]  <- model_sol$parameters$tau/log(2)*
      (log(m0) + (gridMat_reduced/model_sol$parameters$m_pi-m0)/(m0))
    
    lines(gridMat_reduced,linF,
          col=colors[iii],lwd=2,lty=1)
    
    for(j in 1:3){
      if(j == 1){labels <- "RCP4.5"}
      if(j == 2){labels <- "RCP6.0"}
      if(j == 3){labels <- "RCP8.5"}
      rug(matrix_MAT_RCP[j,iii],lwd=2,col=colors[iii],
          ticksize = 0.05*j)
      if(j == 3){
        abline(v=matrix_MAT_RCP[j,iii],lty=3,col=colors[iii])
        text(x=matrix_MAT_RCP[j,iii],y=ylim[1] + 
               ifelse(iii==1,.8,.5)*(ylim[2]-ylim[1]),
             srt=90, pos=2,col=colors[iii],
             labels = paste("RCP8.5, in ",vector.of.dates[iii],sep=""))
      }
      if(j == 2){
        text(x=matrix_MAT_RCP[j,iii],
             y=ylim[1] + 0.05*j*(ylim[2] - ylim[1]),
             srt=0, pos=4, cex=1,col=colors[iii],
             offset = 0.0,
             labels = labels)
      }
    }
  }
  
  legend("topleft",
         legend=c(legend_names,"Nonlinearized"),
         lty=c(rep(1,length(vector.of.dates)),2),
         col=c(colors,"black"),
         lwd=c(2,2),seg.len = 3,
         bty = "n",cex=1.1)
  
}


dev.off()

