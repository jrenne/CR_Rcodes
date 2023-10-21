

multi.lt.fct.Uh_MULTI<-function(model_sol,Uh,X=model_sol$X,t=0){
  H  <-dim(Uh)[3]
  K  <-dim(Uh)[2]
  U  <-list(matrix(Uh[,,H],ncol=K))
  if((t+H)>(model_sol$Tmax-1)){
    a.h<-matrix(rbind(c(a1.fct.inf(model_sol,U[[1]])),
                      matrix(NaN,H-1,K)),H,K)   
    i<-1
    #as long as t+H-i is larger than 98=max(a1/b1), we stay time-indep
    while((t+H-i)>(model_sol$Tmax-2)&i<H){
      U[[i+1]]<-Uh[,,H-i]+b1.fct.inf(model_sol,U[[i]])
      a.h[i+1]<-c(a1.fct.inf(model_sol,U[[i+1]]))
      i<-i+1
    }
    if(i<=(H-1)){
      for (k in i:(H-1)){
        U[[k+1]]<-Uh[,,H-k]+b1.fct(model_sol,U[[k]],t+H-k)
        a.h[k+1]<-c(a1.fct(model_sol,U[[k+1]],t+H-(k+1)))
      } 
    }
  }else{
    a.h<-matrix(rbind(c(a1.fct(model_sol,U[[1]],t+H-1)),
                      matrix(NaN,H-1,K)),H,K)   
    for (k in 1:(H-1)){
      U[[k+1]]<-Uh[,,H-k]+b1.fct(model_sol,U[[k]],t+H-k)
      a.h[k+1,]<- c(a1.fct(model_sol,U[[k+1]],t+H-(k+1)))
    }
  }
  if(t>(model_sol$Tmax-2)){
    uX_t.h<-exp(apply(a.h,2,sum)+t(b1.fct.inf(model_sol,U[[H]]))%*%X)
    psi.0 <-apply(a.h,2,sum)
    psi.1 <-b1.fct.inf(model_sol,U[[H]])
  }else{
    uX_t.h<-exp(apply(a.h,2,sum)+t(b1.fct(model_sol,U[[H]],t))%*%X)
    psi.0 <-apply(a.h,2,sum)
    psi.1 <-b1.fct(model_sol,U[[H]],t) 
  }
  mylist<-list("psi.0"=psi.0,"psi.1"=psi.1,"uX_t.h"=uX_t.h)
  return(mylist)
}

# One layer per maturity:
K <- 3
H <- 10
Uh <- array(0,c(model_sol$n.X,K,H))
Uh[,1,H] <- omega_A
Uh[,2,H] <- 1.1 * omega_A
Uh[,3,] <- .01*rnorm(H*model_sol$n.X)

essai <- multi.lt.fct.Uh_MULTI(model_sol,Uh)

# Compare:
check <- matrix(0,K,1)
for(k in 1:K){
  Uh_single <- matrix(Uh[,k,],model_sol$n.X,H)
  check[k] <- multi.lt.fct.Uh(model_sol,Uh_single)$uX_t.h
}
cbind(essai$uX_t.h,check)



varphi_multi_omega<-function(model_sol,omega,H,X=model_sol$X,t=0){
  # omega of dimension n.X x K
  
  K <- dim(omega)[2]
  
  P.pi <- t(matrix(unlist(model_sol$pi),
                   model_sol$n.X,model_sol$Tmax-2))
  eta0 <- model_sol$eta0
  eta1 <- t(matrix(unlist(model_sol$eta1),
                   model_sol$n.X,model_sol$Tmax-2))
  
  P.a.pi<-matrix(NaN,H-1,1)
  P.b.pi<-matrix(NaN,H-1,model_sol$n.X)
  
  #mu_r0 <- eta0[t+(1:H),]
  #mu_r1 <- eta1[t+(1:H),]
  
  U <- array(0,c(model_sol$n.X,K,H))
  
  for(k in 1:(H-1)){
    P.a.pi[k] <-a1.fct(model_sol,matrix(P.pi[t+k+1,],ncol=1),t+k)
    P.b.pi[k,]<-b1.fct(model_sol,matrix(P.pi[t+k+1,],ncol=1),t+k)
    
    U[,,k] <- - eta1[t+k+1,] - P.b.pi[k,] + P.pi[t+k,]
  }
  U[,,H] <- P.pi[t+H,] + omega
  
  res.psi <- multi.lt.fct.Uh_MULTI(model_sol,Uh=U,X=X,t=t)
  
  # Specific treatment for k=1:
  P.a.pi.1 <- a1.fct(model_sol,matrix(P.pi[t+1,],ncol=1),t)
  P.b.pi.1 <- b1.fct(model_sol,matrix(P.pi[t+1,],ncol=1),t)
  
  varphi0 <- - sum(eta0[(t+1):(t+H)]) - c(P.a.pi.1) - sum(P.a.pi) +
    res.psi$psi.0
  varphi1 <- matrix(- eta1[t+1,] - P.b.pi.1,model_sol$n.X,K) + 
    res.psi$psi.1
  
  Varphi <- exp(varphi0 + t(varphi1)%*%X)
  
  mylist<-list("varphi0"=varphi0,"varphi1"=varphi1,"P.t"=Varphi)

  return(mylist)
}



H <- 10

kk <- 10
omega <- cbind(omega_T.at,omega_ZCB,
               matrix(rnorm(kk*model_sol$n.X),model_sol$n.X,kk))
K <- dim(omega)[2]

tictoc::tic()
res.varphi.multi <- varphi_multi_omega(model_sol,omega,H=H)
tictoc::toc()

tictoc::tic()
check <- matrix(0,K,1)
for(k in 1:K){
  res <- varphi(model_sol,omega.varphi = matrix(omega[,k],ncol=1),
                     H = H)
  check[k] <- res$P.t[H]
}
tictoc::toc()

#cbind(res.varphi.multi$P.t,check)


varphi.hat.fast<-function(model_sol,omega,H,x,a,b,X=model_sol$X,t=0){
  dx<-matrix(x-c(0,x[1:length(x)-1]),length(x),1)
  U <- matrix(omega,model_sol$n.X,length(x))
  U <- U + 1i*matrix(a,ncol=1) %*% matrix(x,nrow=1)
  s1<-varphi_multi_omega(model_sol,U,H,X=X,t=t)$P.t
  fx<-outer(x,b,function(r,c)Im(matrix(s1,ncol=1)*exp(-1i*r*c))/r)*dx
  varphi.hat<-varphi(model_sol,omega,H,X,t)[[3]][H]/2-1/pi*sum(fx)
  return(varphi.hat)
}

varphi.bar.fast<-function(model_sol,omega,H,x,a,b,X=model_sol$X,t=0){
  eps       <-10^-5
  varphi.bar<-(varphi.hat.fast(model_sol,eps*omega,H,x,a,b,X,t)[[1]]-
                 varphi.hat.fast(model_sol,0,H,x,a,b,X,t)[[1]])/eps
  return(varphi.bar)
}


omega <- omega_A
varphi.hat.fast(model_sol,omega,H,x,a,b,X=model_sol$X,t=0)
varphi.bar.fast(model_sol,omega,H,x,a,b,X=model_sol$X,t=0)


