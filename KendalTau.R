library("mvtnorm")
kendtau=function(x){
  n=500
  X=cbind(1, x)
  delta=c(0.1, 0.3)
  #The first bivariate normal component
  beta11=c(0.1,0.3)
  beta12=c(0.1,0.4)
  sigma1=0.5
  I=matrix(c(1,0,0,1), nrow = 2)
  #The second bivariate normal component
  beta21=c(0.3, 0.1)
  beta22=c(0.4, 0.1)
  sigma2=1
  
  pi=cbind(exp(X%*%delta)/(exp(X%*%delta)+1), 1/(exp(X%*%delta)+1))
  

  
  u=matrix(rep(0, n*2), nrow = n)
  for (i in 1:n) {
    #Generate "latent" variable indicating which component
    h=t(apply(pi, 1,  rmultinom, n=1, size=1))
    
    mu1=cbind((X%*%beta11)[,1], (X%*%beta12)[,1])
    mu2=cbind((X%*%beta21)[,1], (X%*%beta22)[,1])
    
    #Sample from the first component
    Z1=rmvnorm(1, mu1,sigma=sigma1*I)
    #sample from the second component
    Z2=rmvnorm(1, mu2,  sigma=sigma2*I)
    
    Z=Z1*h[1]+Z2*h[2]
    
    #Get copula data
    u[i,]=pnorm(Z)
  }
 
  C=function(u){
    pi[1]*(pmvnorm(lower=-Inf, upper=qnorm(u), mean=mu1[1,], sigma=sigma1*I))+
      pi[2]*(pmvnorm(lower=-Inf,upper= qnorm(u),mean=mu2[1,], sigma=sigma2*I))
  }  
  
  tau=4*mean(apply(u, 1, C))-1
  return(tau)
}


xx=seq(0,5, length.out = 1000)
kendtau(0)
Y=sapply(xx, kendtau)
plot(xx,Y, type = 'l')
