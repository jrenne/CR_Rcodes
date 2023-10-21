make_confidence_area <- function(x,y,pdf_xy,p,tol=10^(-2)){
  # Computes the coordinates of a polygon that delineates a confidence
  #     area associated with probability p.
  # The components of x (and of y) must be equally spaced.
  
  h.x <- x[2] - x[1]
  h.y <- y[2] - y[1]
  
  pdf_xy_h2 <- pdf_xy * h.x * h.y # The sum of all entries of pdf_xy_h2 should ne close to 1.
  
  # Use bisection to find limit value of pdf_xy_h2 to get probability p within the area:
  z.up <- max(pdf_xy_h2)
  z.lo <- 0
  error <- 1000
  while(error > tol){
    pdf_xy_h2_aux <- pdf_xy_h2
    pdf_xy_h2_aux[pdf_xy_h2<(z.lo+z.up)/2] <- 0
    aux <- sum(pdf_xy_h2_aux) - p
    if(aux>0){
      z.lo <- (z.lo+z.up)/2
    }else{
      z.up <- (z.lo+z.up)/2
    }
    error <- abs(aux)
    print(error)
  }
  
  M <- pdf_xy_h2_aux
  M[M>0] <- 1 # Then the matrix M is filled only with 0 and 1.
  x.polygon.1 <- NULL
  x.polygon.2 <- NULL
  y.polygon.1 <- NULL
  y.polygon.2 <- NULL
  for(i in 1:length(x)){
    if(sum(M[i,])>0){
      x.polygon.1 <- c(x.polygon.1,x[i])
      x.polygon.2 <- c(x.polygon.2,x[i])
      index.y <- which(M[i,]>0)
      y.polygon.1 <- c(y.polygon.1,y[index.y[1]])
      y.polygon.2 <- c(y.polygon.2,y[index.y[length(index.y)]])
    }
  }
  
  x.polygon <- c(x.polygon.1,rev(x.polygon.2),x.polygon.1[1])
  y.polygon <- c(y.polygon.1,rev(y.polygon.2),y.polygon.1[1])
  
  return(list(
    x.polygon = x.polygon,
    y.polygon = y.polygon
  ))
}

x <- seq(-4,30,by=.01)
y <- seq(-4,4,by=.001)
pdf_xy <- matrix(dnorm(x),ncol=1) %*% matrix(dnorm(y),nrow=1)

p <- .95
res <- make_confidence_area(x,y,pdf_xy,p)

plot(res$x.polygon,res$y.polygon,type="l")


x <- res$x.polygon
y <- res$y.polygon
n <- length(x)
t <- 1:n
ts <- seq(1, n, by = 1/10)
xs <- splinefun(t, x,method = "fmm")(ts)
ys <- splinefun(t, y,method = "fmm")(ts)

plot(x, y,col="white")
polygon(x,y,col="grey", border = NaN)
#lines(xs, ys)


