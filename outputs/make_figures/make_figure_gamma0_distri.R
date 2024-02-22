# ==============================================================================
# Two figures illutrating the gamma-zero distributions
# ==============================================================================

t <- 1:50
w <- .8^t
y <- 0
Y <- y
for(i in t){
  y <- .9*y + w[i]
  Y <- c(Y,y)
}

all.mu <- c(0.0001,.01,.1)


pdf(file="outputs/Figures/Figure_gamma0.pdf",pointsize=11,width=7, height=3)

par(mfrow=c(1,length(all.mu)))
par(plt=c(.1,.95,.18,.88))

all.main <- c(expression(paste("Panel (a) ",mu," = 0.0001",sep="")),
              expression(paste("Panel (b) ",mu," = 0.01",sep="")),
              expression(paste("Panel (c) ",mu," = 0.1",sep="")))

for(iii in 1:length(all.mu)){
  
  mu <- all.mu[iii]
  
  plot(Y,type="l",ylim=c(0,1.4*max(Y)),lwd=2,
       xlab="t+h",ylab="",
       main = all.main[iii])
  
  a  <- 0
  b  <- 1
  
  nb.replic <- 10000
  x <- matrix(0,nb.replic,1)
  X <- x
  for(i in t){
    Z <- rpois(nb.replic,1/mu*(a+b*Y[i+1]))
    x <- rgamma(nb.replic,shape=Z,scale=mu)
    X <- cbind(X,x)
  }
  
  lines(Y)
  
  q <- c(.1,.5,.9)
  Q <- apply(X,2,function(x){quantile(x,q)})
  
  lines(Q[2,],col="red",lty=2,lwd=2)
  
  tt <- 1:(max(t)+1)
  polygon(c(tt,rev(tt)),
          c(Q[1,],rev(Q[3,])),
          col="#18A7B526",border = "dark grey")
  
  if(iii==1){
    legend("topright",
           legend=c(expression(paste(y[t+h]," = ",E[t](x[t+h]),sep="")),
                    expression(paste("Median of ",x[t+h],sep=""))),
           lty=c(1,2),
           col=c("black","red"),
           #pch=c(1,0,3),
           cex=1.2,
           lwd=c(2,2,2),bty = "n")
  }
  
}

dev.off()


#------- FOURIER transform for permafrost-related emissions:
FOURIER<-function(model,x,gamma){
  # works only for i^th variable
  dx <- matrix(x-c(0,x[1:length(x)-1]),length(x),1)
  s1 <- matrix(PSI(model,1i*x),ncol=1)
  fx <- outer(x,gamma,function(r,c)Im(s1[,1]*exp(-1i*r*c))/r)*
    dx[,1]
  f <- 1/2-1/pi*apply(fx,2,sum)
  return(f)
}

PSI <- function(model,u){
  return(exp(model$y * u /(1 - u * model$mu) ))
}


all.h <- c(1,16,36)

pdf(file="outputs/Figures/Figure_gamma0_Damages.pdf",pointsize=11,width=7, height=3)

par(mfrow=c(1,length(all.h)))
par(plt=c(.1,.95,.2,.85))

T <- 2
w <- -.0024+.0037*T

all.proba <- round(100*exp(-all.h*w/model$mu))

gamma <- seq(0.00001,1,length.out=200)
n     <- length(gamma)
reduced.gamma <- .5 * gamma[2:n] + .5 * gamma[1:(n-1)]

x <- exp(seq(-2,2,length.out = 10000)) #grid for Proposition 8 (FOURIER)
x <- seq(.00000001,10000,length.out = 100000)

jjj <- 0
for(h in all.h){
  jjj <- jjj + 1
  
  model <- list(y = h*w,
                mu = .0351)
  
  cdf  <- FOURIER(model,x,gamma)
  cdf <- pmax(pmin(cdf,1),0)
  for(i in 2:length(cdf)){
    if(cdf[i]<cdf[i-1]){
      cdf[i] <- cdf[i-1]
    }
  }

  tmp <- splinefun(x=gamma, y=cdf, method="hyman")
  
  fitted.cdf.values <- tmp(gamma)
  fitted.pdf.values <- diff(fitted.cdf.values) / (gamma[2] - gamma[1])

  plot(reduced.gamma,fitted.pdf.values,type="l",
       ylim=c(0,8),
       xlim=c(0,.6),
       col="white",
       las=1,
       main=paste("Horizon h = ",h," (",5*h," years)",sep=""),
       xlab=expression(paste("Sum of damages ",sum(D[t+i], i==1, h),sep="")))
  
  polygon(c(reduced.gamma,rev(reduced.gamma)),
          c(fitted.pdf.values,rev(0*fitted.pdf.values)),
          col="grey",border="black")
  
  eval(parse(text = gsub(" "," ",
                         paste("Text <- expression(paste('Proba. that ',sum(D[t+i], i==1, h),' = 0 is ',",
                               all.proba[jjj],",'%',
                        sep=''))",sep=""))))
  
  text(.3,5,Text,cex=1)
  text(.3,4.3,
       "(Dirac mass at zero)",
       cex=1)
  
  #print(exp(-model$y/model$mu))
}


dev.off()



