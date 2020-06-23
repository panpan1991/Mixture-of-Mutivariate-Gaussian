library(pscl)
library(MASS)
library(pgdraw)
library(matrixcalc)
library(ggcorrplot)
#load dataset
#dat=read.table("data4hw6")
dat=read.table("dat")
attach(dat)
n=length(Y1)

#construct design matrix
X=matrix(c(rep(1,n),X), ncol = 2)



k=10 #number of components in the mixture of regression
  
#Hyperparameter for sigma's square's prior
alpha=1
lamda=1
#Hyperparameter for deltas' prior
B0=diag(100, nrow=2,ncol=2)

#Starting point
sigma.sq=rep(1,k)
beta=matrix(rep(c(1,1),k),ncol = k)

sigma2.sq=rep(1,k)
beta2=matrix(rep(c(1,1),k),ncol = k)

delta=matrix(rep(c(1,1),k),ncol = k)

#This is an n*k matrix
p=exp(X%*%delta)/apply(exp(X%*%delta), 1, sum)
#initial indicator matrix
Z=t(apply(p, 1, rmultinom, n=1, size=1))


nk=apply(Z, 2, sum)
#Xk[[j]] is partition of design matrix that assign to regression j
Xk=list()
for (j in 1:k) {
  Xk[[j]]=X[Z[,j]==1,]
}

####################Iterations###################################
n.iteration=1000
#containers
sigma.sq.record=matrix(NA, nrow = n.iteration, ncol = k)
beta.record=list()
sigma2.sq.record=matrix(NA, nrow = n.iteration, ncol = k)
beta2.record=list()

delta.record=list()

Z.record=list()
nk.record=list()
Xk.record=list()

 error=0
# i=1
for(i in 1:n.iteration){

#???
individual.density=dnorm(dat$Y1,mean=X%*%beta,
                         sd=matrix(rep(sqrt(sigma.sq),n), byrow = TRUE, ncol = k))


# for (ii in 1:n) {
#   dnorm(Y[ii],mean=X[ii,]%*%beta,  sd=sqrt(sigma.sq))
# }



mixture=apply(p*individual.density, 1, sum)
#mixture=matrix(rep(mixture,3), ncol = 3)
prob=p*individual.density/mixture
Z=t(apply(prob, 1, rmultinom, n=1, size=1))
Z.record[[i]]=Z
###########################################################
nk=apply(Z, 2, sum)
nk.record[[i]]=nk

for (j in 1:k) {
  Xk[[j]]=X[Z[,j]==1,]
}

Xk.record[[i]]=Xk

for (j in 1:k) {
  # Since we have two parameters beta_0 and beta_1 to estimate, 
  # we need two observations being allocated to each expert at least
  if(nk[j]>=2){

    #sigmas
    D=diag(Z[,j])
    sigma.sq[j]=rigamma(1, 
                      sum(Z[,j])/2+alpha+1, 
                      1/2*t(Y1-X%*%beta[,j])%*%D%*%(Y1-X%*%beta[,j])+
                        lamda+1/2*t(beta[,j])%*%((1/nk[j])*t(Xk[[j]])%*%Xk[[j]])%*%beta[,j])
    
    sigma2.sq[j]=rigamma(1, 
                        sum(Z[,j])/2+alpha+1, 
                        1/2*t(Y2-X%*%beta2[,j])%*%D%*%(Y2-X%*%beta2[,j])+
                          lamda+1/2*t(beta2[,j])%*%((1/nk[j])*t(Xk[[j]])%*%Xk[[j]])%*%beta2[,j])
    
    
    #betas
    B=diag(Z[,j]/sigma.sq[j])
    m0=nk[j]*sigma.sq[j]*solve(t(Xk[[j]])%*%Xk[[j]])
    m=solve(t(X)%*%B%*%X+solve(m0))
    beta[,j]=mvrnorm(1,m%*%t(X)%*%B%*%Y1,m)
    
    B=diag(Z[,j]/sigma2.sq[j])
    m0=nk[j]*sigma2.sq[j]*solve(t(Xk[[j]])%*%Xk[[j]])
    m=solve(t(X)%*%B%*%X+solve(m0))
    beta2[,j]=mvrnorm(1,m%*%t(X)%*%B%*%Y2,m)
  }else{
    error=error+1
  }
}

##################################################################
#W=rpg(500,1,X%*%delta)
for (j in 1:k) {
  
  Cj=log(apply(exp(X%*%delta[,-j]), 1, sum))
  
  #W=pgdraw(1, X%*%delta[,j]-Cj)
  W=pgdraw(1, X%*%delta[,j]-Cj)
  K=Z[,j]-1/2
  omega=diag(W)
  mw=solve(t(X)%*%omega%*%X+solve(B0))%*%(t(X)%*%(K+omega%*%Cj))
  vw=solve(t(X)%*%omega%*%X+solve(B0))
  
  delta[,j]=mvrnorm(1,mw,vw)
}
p=exp(X%*%delta)/apply(exp(X%*%delta), 1, sum)

sigma.sq.record[i,]=sigma.sq
beta.record[[i]]=beta
sigma2.sq.record[i,]=sigma.sq
beta2.record[[i]]=beta
delta.record[[i]]=delta

print( i)

}


#Trace plots
par(mfrow=c(3,2))
beta1.record=matrix(rep(0,2*n.iteration), nrow = n.iteration)
for (i in 1:n.iteration) {
  beta1.record[i,]=beta.record[[i]][,1]
}

beta2.record=matrix(rep(0,2*n.iteration), nrow = n.iteration)
for (i in 1:n.iteration) {
  beta2.record[i,]=beta.record[[i]][,2]
}

beta3.record=matrix(rep(0,2*n.iteration), nrow = n.iteration)
for (i in 1:n.iteration) {
  beta3.record[i,]=beta.record[[i]][,3]
}


plot(beta1.record[,2], main="Beta1 slope")
plot(beta1.record[,1], main="Beta1 intercept")

plot(beta2.record[,2], main="Beta2 slope")
plot(beta2.record[,1], main="Beta2 intercept")

plot(beta3.record[,2], main="Beta2 slope")
plot(beta3.record[,1], main="Beta2 intercept")

delta1.record=matrix(rep(0,2*n.iteration), nrow = n.iteration)
for (i in 1:n.iteration) {
  delta1.record[i,]=delta.record[[i]][,1]
}

plot(delta1.record[,2], main="delta1 slope")
plot(delta1.record[,1], main="delta1 intercept")

delta1.record=matrix(rep(0,2*n.iteration), nrow = n.iteration)
for (i in 1:n.iteration) {
  delta1.record[i,]=delta.record[[i]][,1]
}


par(mfrow=c(1,3))
plot(sigma.sq.record[,1], main="Sigma 1")
plot(sigma.sq.record[,2], main="Sigma 2")
plot(sigma.sq.record[,3], main="Sigma 2")

par(mfrow=c(2,1))

#Fitted Curve
p.record=list()
for (i in 1:n.iteration) {
  p.record[[i]]=exp(X%*%delta.record[[i]])/apply(exp(X%*%delta.record[[i]]), 1, sum)
}
fitted=list()
for (i in 1:n.iteration) {
  fitted[[i]]=apply(p.record[[i]]*(X%*%beta.record[[i]]), 1,sum)
}


summ=rep(0,n)
for (i in 500:n.iteration) {
  summ=summ+fitted[[i]]
}

fit=summ/(n.iteration-500)

plot(X[,2],Y1, xlab="X", ylab = "Y1", main = "Sample points and fitted line")
lines(X[,2], fit, col="red", type = "l",lwd=3)
#lines(X[,2],h3(X[,2]), col="green", type = 'l', lwd=3)

# Add a legend
legend(0.8, 3, legend=c("fitted", "h3"),
       col=c("red", "green"), lty=1:1, lwd=3:3, cex=1.5)


#######################

#Fitted Curve
p.record=list()
for (i in 1:n.iteration) {
  p.record[[i]]=exp(X%*%delta.record[[i]])/apply(exp(X%*%delta.record[[i]]), 1, sum)
}
fitted=list()
for (i in 1:n.iteration) {
  fitted[[i]]=apply(p.record[[i]]*(X%*%beta2.record[[i]]), 1,sum)
}


summ=rep(0,n)
for (i in 500:n.iteration) {
  summ=summ+fitted[[i]]
}

fit=summ/(n.iteration-500)

plot(X[,2],Y2, xlab="X", ylab = "Y2", main = "Sample points and fitted line")
lines(X[,2], fit, col="red", type = "l",lwd=3)
#lines(X[,2],h3(X[,2]), col="green", type = 'l', lwd=3)

# Add a legend
legend(0.8, 3, legend=c("fitted", "h3"),
       col=c("red", "green"), lty=1:1, lwd=3:3, cex=1.5)




sigma.sq.record[1000,]
individual.mu=X%*%beta.record[[1000]]
individual.sigma=matrix(rep(sigma.sq.record[1000,], 1500), byrow = TRUE, nrow = 1500)
individual.pred=apply(individual.mu, 1, rnorm, n=10,sd=sigma.sq.record[1000,])

h=t(apply(p.record[[1000]], 1,  rmultinom, n=1, size=1))

Y1.pred=apply(h*t(individual.pred), 1, sum)

plot(X[,2],Y1)
par(new=TRUE)
points(X[,2],Y1.pred, col='green')
###################################
sigma2.sq.record[1000,]
individual.mu2=X%*%beta2.record[[1000]]
individual.sigma2=matrix(rep(sigma2.sq.record[1000,], 1500), byrow = TRUE, nrow = 1500)
individual.pred2=apply(individual.mu2, 1, rnorm, n=10,sd=sigma2.sq.record[1000,])

#h=t(apply(p.record[[1000]], 1,  rmultinom, n=1, size=1))

Y2.pred=apply(h*t(individual.pred2), 1, sum)

plot(X[,2],Y2)
par(new=TRUE)
points(X[,2],Y2.pred, col='green')

Z.pred=cbind(Y1.pred, Y2.pred)
u.pred=t(apply(Z, 1, pnorm))
U.pred=u[,1]
V.pred=u[,2]

par(mfrow=c(1,3))
plot(U[1:(1/3*n)], V[1:(1/3*n)], xlab = 'U', ylab='V', main = '0<z<-1.5')
points(U.pred[1:(1/3*n)], V.pred[1:(1/3*n)], col='red')
plot(U[(1/3*n):(2/3*n)], V[(1/3*n):(2/3*n)],xlab = 'U', ylab='V', main = '1.5<z<3')
points(U.pred[(1/3*n):(2/3*n)], V.pred[(1/3*n):(2/3*n)], col='red')
plot(U[(2/3*n):n], V[(2/3*n):n], xlab = 'U', ylab='V', main = '3<z<5')
points(U.pred[(2/3*n):n], V.pred[(2/3*n):n], col='red')









