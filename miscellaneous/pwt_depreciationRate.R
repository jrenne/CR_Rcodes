rm(list=ls(all=T)) 

library("pwt10")

data("pwt10.01")

delta <- mean(pwt10.01$delta,na.rm=TRUE)*100
densityDelta <- density(pwt10.01$delta,na.rm=TRUE)
plot(densityDelta$x,densityDelta$y,type="l")
var <- var(pwt10.01$delta,na.rm=TRUE)*100

bounds.up <- delta + 2 * sqrt(var)