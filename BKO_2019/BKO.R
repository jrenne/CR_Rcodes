#clear environment
rm(list=ls(all=T)) 

source("procedures/functions_other_models.R")

model.BKO <- make.BKO.model()
model_sol.BKO <- solve_model.bko(model.BKO)
Tbar <- model.BKO$Tbar

H <- 80

res_prices <- compute_ab.bko(model_sol.BKO,H)

plot(res_prices$all_r_a + res_prices$all_r_b*Tbar,type="l")

nb.replic <- 10
res.simul <- simul.model.BKO(model.BKO,
                             H=H,nb.replic=nb.replic)
par(mfrow=c(2,1))
variable <- res.simul$all.delc
ylim <- c(min(variable),max(variable))
plot(1:H,1:H,ylim=ylim,col="white")
for(i in 1:nb.replic){
  lines(variable[,i],type="l",col=i)
}
variable <- res.simul$all.D
ylim <- c(min(variable),max(variable))
plot(1:H,1:H,ylim=ylim,col="white")
for(i in 1:nb.replic){
  lines(variable[,i],type="l",col=i)
}

# Price of Temperature swap:
epsilon <- .0001

B <- exp(res_prices$all_a + res_prices$all_b*Tbar)
res_prices_eps <- compute_ab.bko(model_sol.BKO,H,u=epsilon)

Tswap <- ((res_prices_eps$all_a - res_prices$all_a) +
            (res_prices_eps$all_b - res_prices$all_b)*Tbar)/epsilon
plot(rep(Tbar,H),ylim=c(min(Tbar,Tswap),max(Tbar,Tswap)),type="l")
lines(Tswap,col="red")


