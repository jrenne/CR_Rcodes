# ==============================================================================
# Response of temperature conditional on some MAT paths
# ==============================================================================


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

# # Check ability to replicate ACE results on RCPs:
# RF <- RCP[,14:17]
# RF <- cbind(RF_RCP45,RF_RCP60,RF_RCP85)
# Temp.ini <- matrix(c(model_sol$vector.ini$ini_Tat,
#                      model_sol$vector.ini$ini_Tlo,
#                      0),3,dim(RF)[2])
# forcing <- t(exp(log(2) / eta * RF))
# res <- TempSimulation_ACE(Temp.ini, forcing)
# temp_ace <- res[1,,]
# plot(years,temp_ace[2,])
# lines(RCP_ACE$V1,T_ACE_RCP60)

# Compute T_AT with ACE under CR -----------------------------------------------
# Compute E(M_AT) under CR baseline model:
EV <- EV.fct(model_sol)
RCP_CR <- c(model_sol$vector.ini$ini_Mat,
            EV$EX$M_at)
years_CR <- c(2020,EV$date)
indicators_CR <- which(years_CR %in% years) # to find appropriate years
RCP_CR <- RCP_CR[indicators_CR]
# Compute weight in RCP45 (versus RCP60) to approximate MAT in CR with
# linear combination of RCP45 and RCP60:
coeff.45  = (mean(RCP_CR) - mean(MAT_RCP60))/(mean(MAT_RCP45)-mean(MAT_RCP60))
RF_CR_4ACE <- coeff.45 * RF_RCP45 + (1 - coeff.45) * RF_RCP60 # to be used in TempSimulation_ACE
Temp.ini <- matrix(c(model_sol$vector.ini$ini_Tat,
                     model_sol$vector.ini$ini_Tlo,0),3,2)
eta <- 3.8
forcing <- matrix(1,2,1) %*% t(exp(log(2) / eta * RF_CR_4ACE))
res <- TempSimulation_ACE(Temp.ini, forcing)
T_ACE_CR <- res[1,1,]
# ------------------------------------------------------------------------------






FILE = "/outputs/Figures/Figure_RCP_to_TAT.pdf"
pdf(file=paste(getwd(),FILE,sep=""),pointsize=12, width=7, height=7)

#par(mfrow=c(2,2))
par(plt=c(.06,.97,.16,.86))

nf <- layout(
  matrix(c(1,1,1,2,3,4), ncol=3, byrow=TRUE), 
  widths=c(1,1,1), 
  heights=c(3,2)
)


# Grid of values of M_AT (for first and second panel):
gridMat<-seq(900,3200,by=1)

nlinF <- matrix(0,length(gridMat),1)
nlinF[1:length(nlinF)] <- model_sol$parameters$tau/log(2)*
  log(gridMat[1:length(gridMat)]/model_sol$parameters$m_pi)

vector.of.dates <- c(2070,2100,2150)
colors <- c("#CCCCCC","#999999","#666666")
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
     main=expression(paste("(a) Relationship between radiative forcings and atmospheric carbon concentration")))
grid()

title(ylab=expression(paste("FCO"[2]," (in Wm-2)")), line=1.5,
      cex.lab=1)


# Distributions of M_AT for CR model -------------------------------------------
# For Fourier-transform computation:
x <- exp(seq(-10,10,length.out = 5000)) # grid for Proposition 8 (Fourier)
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
  
  #x=M_at, y=F
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
  #print(quantile_RCP85_considered)
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
    if(j == 2){labels <- "RCP6"}
    if(j == 3){labels <- "RCP8.5"}
    rug(matrix_MAT_RCP[j,iii],lwd=2,col=colors[iii],
        ticksize = 0.05*j)
    if(iii==2){
      text(x=matrix_MAT_RCP[j,iii],
           y=ylim[1] + 0.05*j*(ylim[2] - ylim[1]),
           srt=0, pos=4, cex=.8,
           offset = 0.0,
           labels = labels)
    }
  }
}

legend("topleft",
       legend=c(legend_names,"Non-linearized"),
       lty=c(rep(1,length(vector.of.dates)),2),
       col=c(colors,"black"),
       lwd=c(2,2),seg.len = 3,
       bty = "n",cex=1.1)




#determine  quantiles of RCP8.5 values (for analysis):


# # RCP 4.5:
# abline(v=MAT_RCP45[indic.first.in.yearsvector],col="black",
#        lty=1,lwd=1)
# text(x=MAT_RCP45[indic.first.in.yearsvector],y=.98*ymax,
#      srt=90, pos=2,
#      labels = paste("RCP4.5, in 2100 and ",other.year,sep=""))
# abline(v=MAT_RCP45[indic.other.in.yearsvector],col="black",
#        lty=1,lwd=2)
# 
# # RCP 6.5:
# abline(v=MAT_RCP60[indic.first.in.yearsvector],col="black",
#        lty=1,lwd=1)
# text(x=MAT_RCP60[indic.first.in.yearsvector],y=.98*ymax,
#      srt=90, pos=2,
#      labels = "RCP6, in 2100")
# abline(v=MAT_RCP60[indic.other.in.yearsvector],col="black",
#        lty=1,lwd=2)
# text(x=MAT_RCP60[indic.other.in.yearsvector],y=.98*ymax,
#      srt=90, pos=2,
#      labels = paste("RCP6, in ",other.year,sep=""))
# 
# # RCP 8.5:
# abline(v=MAT_RCP85[indic.first.in.yearsvector],col="black",
#        lty=1,lwd=1)
# text(x=MAT_RCP85[indic.first.in.yearsvector],y=.98*ymax,
#      srt=90, pos=2,
#      labels = "RCP8.5, in 2100")
# abline(v=MAT_RCP85[indic.other.in.yearsvector],col="black",
#        lty=1,lwd=2)
# text(x=MAT_RCP85[indic.other.in.yearsvector],y=.98*ymax,
#      srt=90, pos=2,
#      labels = paste("RCP8.5, in ",other.year,sep=""))

# MAT4arrow1 <- 1550
# MAT4arrow2 <- 1900
# y4arrow1 <- Mat.pdf.first[which(gridMat[-1]==MAT4arrow1)]
# y4arrow2 <- Mat.pdf.other[which(gridMat[-1]==MAT4arrow2)]
# 
# arrowlength <- 800
# upshift <- .2*ymax
# 
# arrows(x0 = MAT4arrow1 + arrowlength, y0 = y4arrow1 + upshift,
#        x1 = MAT4arrow1, y1 = y4arrow1,angle = 20,
#        length = 0.1, col = "dark grey", lwd = 1)
# text(x = MAT4arrow1 + arrowlength, y = y4arrow1 + upshift,
#      pos=4,
#      labels = expression(paste("Model-implied p.d.f. of ",M[AT]," in 2100",sep="")))
# 
# arrows(x0 = MAT4arrow2 + arrowlength, y0 = y4arrow2 + upshift,
#        x1 = MAT4arrow2, y1 = y4arrow2,angle = 20,
#        length = 0.1, col = "dark grey", lwd = 2)
# text(x = MAT4arrow2 + arrowlength, y = y4arrow2 + upshift,
#      pos=4,
#      labels = expression(paste("Model-implied p.d.f. of ",M[AT]," in 2200",sep="")))





# Last plots -------------------------------------------------------------------


par(plt=c(.25,.95,.15,.85))

for(rcp in c(45,60,85)){
  if(rcp==00){
    Mat.trajectory <- RCP_CR[-1]
    Tat.ACE        <- T_ACE_CR
    main.t =expression(paste("(XX) ",E(M[AT])," path from present framework",sep=""))
  }
  if(rcp==45){
    Mat.trajectory <- MAT_RCP45[-1]
    Tat.MAGICC     <- T_RCP45
    Tat.ACE        <- T_ACE_RCP45
    main.t =expression(paste("(b) ",T[AT]," resulting from RCP4.5",sep=""))
  }
  if(rcp==60){
    Mat.trajectory <- MAT_RCP60[-1]
    Tat.MAGICC     <- T_RCP60
    Tat.ACE        <- T_ACE_RCP60
    main.t =expression(paste("(c) ",T[AT]," resulting from RCP6",sep=""))
  }
  if(rcp==85){
    Mat.trajectory <- MAT_RCP85[-1]
    Tat.MAGICC     <- T_RCP85
    Tat.ACE        <- T_ACE_RCP85
    main.t =expression(paste("(d) ",T[AT]," resulting from RCP8.5",sep=""))
  }
  
  res_CR    <- simul_TAT_condit_MAT(model_sol,Mat.trajectory)
  res_CDICE <- simul_TAT_condit_MAT_CDICE(model_sol$tstep,Mat.trajectory,
                                          model_sol$vector.ini$ini_Tat,
                                          model_sol$vector.ini$ini_Tlo,
                                          model_sol$vector.ini$ini_F,
                                          model_sol$parameters)
  Tat.nonlinear <- res_CR$Tat.nonlinear
  Tat.linear    <- res_CR$Tat.linear
  
  plot(years,Tat.nonlinear,type="l",lwd=2,col="white",
       ylim=c(1,1.25*max(Tat.nonlinear,Tat.linear)),
       las=1,xlab="",ylab="",
       main=main.t)
  grid()
  
  title(ylab='Atm Temperature Anomaly, in °C', line=3,
        cex.lab=1)
  
  lines(years,Tat.nonlinear,col="dark grey",lwd=5,pch=3)
  lines(years,Tat.linear,col="black",lwd=4,lty=3)
  lines(years,res_CDICE$Tat,col="red",lwd=1)
  
  if(rcp==00){
    lines(years,Tat.ACE,col="#006400",lty=4,lwd=2)
  }else{
    lines(RCP_ACE$V1,Tat.ACE,col="#006400",lty=4,lwd=2)
  }
  
  if(rcp==45){
    legend("topleft",
           legend=c("CR Linearized",
                    "CR Non-linearized",
                    "CDICE",
                    "ACE"),
           lty=c(3,1,1,4),
           col=c("black","dark grey","red","#006400"),
           pch=c(NaN),
           lwd=c(4,4,1,2),
           seg.len = 4,
           bty = "n",cex=1)
  }
}


dev.off()

