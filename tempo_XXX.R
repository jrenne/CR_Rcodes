

RCP_MAGICC <- read.csv("data/RCP_Mat_MAGICC.csv", header=FALSE)

RF <- RCP_MAGICC[,14:17]
indicators <- seq(255,636,by=5)
RF <- RF[indicators,]

Temp.ini <- rbind(
  c(1.0231,    1.0122,    0.9890,    1.0265),
  c(0.2781,    0.2780,    0.2768,    0.2787),
  c(0.0131,    0.0131,    0.0130,    0.0131))

forcing <- t(exp(log(2) / eta * RF))

res <- TempSimulation_ACE(Temp.ini, forcing)

TempSim <- res[1,,]
plot(TempSim[3,])

# see paper (A.18)
sigma5yr = matrix(c(0.0000104,
                    0.036187,
                    0,0.486022,0.953354,
                    0.001918,0,0.01045,0.998081),3,3)
sigma_forc <- 0.513966
Mpre <- 596.4 # from Matlab codes

# Determine Gt process in Traeger (2023), using his eq.(6):
eta <- 3.8
Mpre <- 596.4 # from Matlab codes
RF <- RCP_MAGICC[,14:17]
M1 <- RCP_MAGICC[,2:5]
G <- exp(RF*log(2)/eta)*Mpre - M1
plot(G[,4])
lines(G[,3])

# Determine E(M1) in CR:
EV <- EV.fct(model_sol)
EM1_CR <- EV$EX$M_at
plot(EV$date,EM1_CR,ylim=c(0,2000))
indicators_in_CR <- which(RCP_MAGICC$V1 %in% EV$date)
M1_RCP45 <- M1[indicators_in_CR,2]
M1_RCP60 <- M1[indicators_in_CR,3]
lines(EV$date,M1_RCP45)
lines(EV$date,M1_RCP60)
coeff.45  = (mean(EM1_CR) - mean(M1_RCP60))/(mean(M1_RCP45)-mean(M1_RCP60))
lines(EV$date,coeff.45 * M1_RCP45 + (1 - coeff.45) * M1_RCP60,col="red")

  
