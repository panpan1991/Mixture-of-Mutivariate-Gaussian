library("ggplot2")
library("mixtools")
library("mvtnorm")

set.seed(100)
n=1000
X=seq(0,5, length.out = n)
X=cbind(rep(1, n), X)
delta=c(0, 1)

pi=cbind(exp(X%*%delta)/(exp(X%*%delta)+1), 1/(exp(X%*%delta)+1))

plot(X[,2], pi[,1])
plot(X[,2], pi[,2])

#Generate "latent" variable indicating which component
h=t(apply(pi, 1,  rmultinom, n=1, size=1))


#The first bivariate normal component
beta11=c(1,0.5)
beta12=c(1,2)
sigma1=0.5
I=matrix(c(1,0,0,1), nrow = 2)
mu1=cbind((X%*%beta11)[,1], (X%*%beta12)[,1])

plot(X[,2], mu1[,1])
plot(X[,2], mu1[,2])

#Sample from the first component
Z1=t(apply(mu1, 1, rmvnorm, n=1, sigma=sigma1*I))



#The second bivariate normal component
beta21=c(1, -5)
beta22=c(1, -0.3)
sigma2=5
I=matrix(c(1,0,0,1), nrow = 2)
#mean vector
mu2=cbind((X%*%beta21)[,1], (X%*%beta22)[,1])
plot(X[,2], mu2[,1])
plot(X[,2], mu2[,2])
#sample from the second component
Z2=t(apply(mu2, 1, rmvnorm, n=1, sigma=sigma2*I))

#use h to select from Z1 and Z2
Z=matrix(rep(0, 2*n), nrow = n)
for (i in 1:n) {
Z[i,]=Z1[i,]*h[i,1]+Z2[i,]*h[i,2]
}

plot(Z[,1], Z[,2])
plot(X[,2], Z[,1])
plot(X[,2], Z[,2])

dat=data.frame(Z[,1], Z[,2], X[,2])
colnames(dat)=c("Y1", "Y2", "X")
write.table(dat, "dat")

#Get copula data
u=t(apply(Z, 1, pnorm))
U=u[,1]
V=u[,2]

par(mfrow=c(1,3))
plot(U[1:(1/3*n)], V[1:(1/3*n)], xlab = 'U', ylab='V', main = '0<z<-1.5')
plot(U[(1/3*n):(2/3*n)], V[(1/3*n):(2/3*n)],xlab = 'U', ylab='V', main = '1.5<z<3')
plot(U[(2/3*n):n], V[(2/3*n):n], xlab = 'U', ylab='V', main = '3<z<5')


#Set up marginal distribution
R=6
mu=c(-9, -5.4, -1.8, 1.8, 5.4, 9)
sigmar=rep(1/sqrt(10), R)
pi=rep(1/R, R)

# evaluate the function at the point x, where the components 
# of the mixture have weights w, means stored in u, and std deviations
# stored in s - all must have the same length.
F = function(x,w,u,s) sum( w*pnorm(x,mean=u,sd=s) )

#Marginal quantile function
# provide an initial bracket for the quantile. default is c(-1000,1000). 
F_inv = function(p,w,u,s,br=c(-1000,1000))
{
  G = function(x) F(x,w,u,s) - p
  return( uniroot(G,br)$root ) 
}


#Get the orginal data
X1=sapply(V, F_inv, w=pi, u=mu, s=sigmar)
X2=sapply(U, F_inv, w=pi, u=mu, s=sigmar)

par(mfrow=c(1,3))
plot(X1[1:(1/3*n)], X2[1:(1/3*n)],  xlab = 'X', ylab='Y', main = '0<z<-1.5')
plot(X1[(1/3*n):(2/3*n)], X2[(1/3*n):(2/3*n)], xlab = 'X', ylab='Y', main = '1.5<z<3')
plot(X1[(2/3*n):n], X2[(2/3*n):n], xlab = 'X', ylab='Y', main = '3<z<5')



######################################################
#Conditional Clayton copula

n=500
U1=runif(n)
V=runif(n)

X=seq(2,5, length.out = n)
kendallTau=rep(0,n)
U2=rep(0,n)
for (i in 1:n) {
  theta=exp(0.8*X[i] - 2)
  kendallTau[i]=theta/(theta+2)
  
  U2[i]=(U1[i]^(-theta)*(V[i]^(-theta/(1+theta))+1))^(-1/theta)
}

par(mfrow=c(1,3))
plot(U1[1:(1/3*n)], U2[1:(1/3*n)], xlab = 'U1', ylab='U2', main = '2<x<-3')
plot(U1[(1/3*n):(2/3*n)], U2[(1/3*n):(2/3*n)],xlab = 'U1', ylab='U2', main = '3<x<4')
plot(U1[(2/3*n):n], U2[(2/3*n):n], xlab = 'U1', ylab='U2', main = '4<x<5')
par(mfrow=c(1,1))
plot(X, kendallTau, type='l')

dat=data.frame(qnorm(U1), qnorm(U2), X)
colnames(dat)=c("Y1", "Y2", "X")
write.table(dat, "dat")
plot(dat$X, dat$Y1)
plot(dat$X, dat$Y2)
#Conditional Frank copula

n=500
U1=runif(n)
V=runif(n)

X=seq(2,5, length.out = n)

fun=function(t){
  t/(exp(t)-1)
}

D=function(theta){
  1/theta*integrate(fun, 0, theta)$value
}
kendallTau=rep(0,n)
U2=rep(0,n)
for (i in 1:n) {
  theta=exp(0.8*X[i] - 2)
  kendallTau[i]=1+4/theta*(D(theta)-1)
  
  U2[i]=-1/theta*log(
    1+(V[i]*(1-exp(-theta)))/
      (V[i]*(exp(-theta*U1[i])-1)-exp(-theta*U1[i]))
                     )
}

par(mfrow=c(1,3))
plot(U1[1:(1/3*n)], U2[1:(1/3*n)], xlab = 'U1', ylab='U2', main = '2<x<-3')
plot(U1[(1/3*n):(2/3*n)], U2[(1/3*n):(2/3*n)],xlab = 'U1', ylab='U2', main = '3<x<4')
plot(U1[(2/3*n):n], U2[(2/3*n):n], xlab = 'U1', ylab='U2', main = '4<x<5')
par(mfrow=c(1,1))
plot(X, kendallTau, type='l')
######################
#Conditional Gumbel copula
n=500
V1=runif(n)
V2=runif(n)

K<-function(w){
  w*(1-log(w)/theta)-V2[i]
}
X=seq(2,5, length.out = n)

U2=rep(0,n)
for (i in 1:n) {
  theta=exp(0.8*X[i] - 2)
  w=uniroot(K, lower = 0, upper = 1)$root
  U1[i]=exp(V1^(1/theta)*log(w))
  U2[i]=exp((1-V1)^(1/theta)*log(w))
}

par(mfrow=c(1,3))
plot(U1[1:(1/3*n)], U2[1:(1/3*n)], xlab = 'U1', ylab='U2', main = '2<x<-3')
plot(U1[(1/3*n):(2/3*n)], U2[(1/3*n):(2/3*n)],xlab = 'U1', ylab='U2', main = '3<x<4')
plot(U1[(2/3*n):n], U2[(2/3*n):n], xlab = 'U1', ylab='U2', main = '4<x<5')

################

Y1=qnorm(U)
Y2=qnorm(V)
dat=data.frame(Y1, Y2, X)
colnames(dat)=c("Y1", "Y2", "X")
write.table(dat, "dat")

#Set up marginal distribution
R=6
mu=c(-9, -5.4, -1.8, 1.8, 5.4, 9)
sigmar=rep(1/sqrt(10), R)
pi=rep(1/R, R)

# evaluate the function at the point x, where the components 
# of the mixture have weights w, means stored in u, and std deviations
# stored in s - all must have the same length.
F = function(x,w,u,s) sum( w*pnorm(x,mean=u,sd=s) )

#Marginal quantile function
# provide an initial bracket for the quantile. default is c(-1000,1000). 
F_inv = function(p,w,u,s,br=c(-1000,1000))
{
  G = function(x) F(x,w,u,s) - p
  return( uniroot(G,br)$root ) 
}

X1=sapply(V, F_inv, w=pi, u=mu, s=sigmar)
X2=sapply(U, F_inv, w=pi, u=mu, s=sigmar)

par(mfrow=c(1,3))
plot(X1[1:(1/3*n)], X2[1:(1/3*n)],  xlab = 'X', ylab='Y', main ='-4<z<-1.3')
plot(X1[(1/3*n):(2/3*n)], X2[(1/3*n):(2/3*n)], xlab = 'X', ylab='Y', main =  '-1.3<z<1.3')
plot(X1[(2/3*n):n], X2[(2/3*n):n], xlab = 'X', ylab='Y', main = '1.3<z<4')


copula=data.frame(V, U)
# write.csv(data, "UV.csv", row.names = FALSE)

Z=qnorm(U);
W=qnorm(V);
data=data.frame(Z, W)
# plot(Z,W)
# write.csv(data, "ZW.csv", row.names = FALSE)

########################################################################


#########################################################################
# EM algorithm to estimate a mixture of multivariate normal distributions
est=mvnormalmixEM(data, lambda = NULL, mu = NULL, sigma = NULL, k = 3,
              arbmean = TRUE, arbvar = TRUE, epsilon = 1e-08, 
              maxit = 10000, verb = FALSE)

#########################################################################
#Sampling from the posterior predictive distribution
x.1<-rmvnorm(round(est$lambda[1],2)*1000, est$mu[[1]], est$sigma[[1]])
x.2 <- rmvnorm(round(est$lambda[2],2)*1000, est$mu[[2]], est$sigma[[2]])
x.3 <-rmvnorm(round(est$lambda[3],2)*1000, est$mu[[3]], est$sigma[[3]])

X.1 <- rbind(x.1, x.2, x.3)
pred=data.frame(X.1)
colnames(pred)=c("Z", "W")

copula.pred=data.frame(pnorm(pred$W), pnorm(pred$Z))
colnames(copula.pred)=c("V", "U")

##########################################################################
ggplot(copula) +
  geom_point(aes(x=V,y=U, color='blue'))+
  geom_point(data=copula.pred,aes(x=V,y=U, color = 'red')) +
  scale_color_discrete(name='', labels=c('Pred','Raw data'))

par(mfrow=c(1,1))
plot(data)
par(new=TRUE)
points(pred, col='green')

#Plot of raw data and sampling from predictive distribution
ggplot(data) +
  geom_point(aes(x=Z,y=W, color='blue'))+
  geom_point(data=pred,aes(x=Z,y=W, color = 'red')) +
  scale_color_discrete(name='', labels=c('Pred','Raw'))


####################################

Z=rmvnorm(500, c(0,0), matrix(c(1, 0, 0, 1), nrow = 2))
plot(Z[, 1], Z[, 2])
plot(pnorm(Z[, 1]), pnorm(Z[, 2]))

