

aD <- model_sol$parameters$a_D
bD <- model_sol$parameters$b_D
muD <- model_sol$parameters$mu_D

TAT <- seq(0,4,by=.1)
plot(TAT,1-exp(-aD/muD - bD/muD*TAT))

1-exp(-aD/muD - bD/muD*c(1,2,3,4))

(aD + bD*c(1,2,3,4))/(1-exp(-aD/muD - bD/muD*c(1,2,3,4)))





