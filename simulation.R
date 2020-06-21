library("ggplot2")
library("mixtools")
library("mvtnorm")

n=1500
X=seq(0,5, length.out = n)
X=cbind(rep(1, n), X)
beta1=c(0.1, 0.3)

pi=cbind(exp(X%*%beta1)/(exp(X%*%beta1)+1), 1/(exp(X%*%beta1)+1))

h=t(apply(pi, 1,  rmultinom, n=1, size=1))

Z1=matrix(rep(0, 2*n), nrow = n)

lambda11=c(0.1,0.3)
lambda12=c(0.1,0.4)
sigma=0.5
I=matrix(c(1,0,0,1), nrow = 2)
mu1=cbind((X%*%lambda11)[,1], (X%*%lambda12)[,1])

for (i in 1:n) {
  Z1[i, ]=rmvnorm(1,mu1[i,], sigma*I)
}


Z2=matrix(rep(0, 2*n), nrow = n)

lambda21=c(0.3, 0.1)
lambda22=c(0.4, 0.1)
sigma2=1
I=matrix(c(1,0,0,1), nrow = 2)
mu2=cbind((X%*%lambda21)[,1], (X%*%lambda22)[,1])

for (i in 1:n) {
  Z2[i, ]=rmvnorm(1,mu2[i,], sigma2*I)
}
Z=matrix(rep(0, 2*n), nrow = n)
for (i in 1:n) {
Z[i,]=Z1[i,]*h[i,1]+Z2[i,]*h[i,2]
}

u=t(apply(Z, 1, pnorm))
U=u[,1]
V=u[,2]

par(mfrow=c(1,3))
plot(U[1:(1/3*n)], V[1:(1/3*n)], xlab = 'U', ylab='V', main = '0<z<-1.5')
plot(U[(1/3*n):(2/3*n)], V[(1/3*n):(2/3*n)],xlab = 'U', ylab='V', main = '1.5<z<3')
plot(U[(2/3*n):n], V[(2/3*n):n], xlab = 'U', ylab='V', main = '3<z<5')



X1=sapply(V, F_inv, w=pi, u=mu, s=sigmar)
X2=sapply(U, F_inv, w=pi, u=mu, s=sigmar)

par(mfrow=c(1,3))
plot(X1[1:(1/3*n)], X2[1:(1/3*n)],  xlab = 'X', ylab='Y', main = '0<z<-1.5')
plot(X1[(1/3*n):(2/3*n)], X2[(1/3*n):(2/3*n)], xlab = 'X', ylab='Y', main = '1.5<z<3')
plot(X1[(2/3*n):n], X2[(2/3*n):n], xlab = 'X', ylab='Y', main = '3<z<5')


######################################################
set.seed(100)
R=6
mu=c(-9, -5.4, -1.8, 1.8, 5.4, 9)
sigmar=rep(1/sqrt(10), R)
pi=rep(1/R, R)
theta=40

# evaluate the function at the point x, where the components 
# of the mixture have weights w, means stored in u, and std deviations
# stored in s - all must have the same length.
F = function(x,w,u,s) sum( w*pnorm(x,mean=u,sd=s) )

# provide an initial bracket for the quantile. default is c(-1000,1000). 
F_inv = function(p,w,u,s,br=c(-1000,1000))
{
  G = function(x) F(x,w,u,s) - p
  return( uniroot(G,br)$root ) 
}


l=5
#function alpha(x)
a<-function(x){
  2/3-1/4*(exp(l*x)/(exp(l*x)+1))
}
#function beta(x)
b<-function(x){
  1/2+1/4*(exp(l*x)/(exp(l*x)+1))
}

n=1500
V=runif(n)
W=runif(n)

X=seq(-4,4, length.out = n)

C <- function(u) {
  (1-beta)*u^(1-alpha)*v^(-beta)*(u^(-theta*alpha)+v^(-theta*beta)-1)^(-1/theta)+
    beta*u^(1-alpha)*v^(-beta*(1+theta))*(u^(-theta*alpha)+v^(-theta*beta)-1)^(-1/theta-1)-w
}

U=rep(0,n)
for (i in 1:n) {
  v=V[i]
  w=W[i]
  alpha=a(X[i])
  beta=b(X[i])
  U[i]=uniroot(C, lower = 0, upper = 1)$root
}

par(mfrow=c(1,3))
plot(U[1:(1/3*n)], V[1:(1/3*n)], xlab = 'U', ylab='V', main = '-4<z<-1.3')
plot(U[(1/3*n):(2/3*n)], V[(1/3*n):(2/3*n)],xlab = 'U', ylab='V', main = '-1.3<z<1.3')
plot(U[(2/3*n):n], V[(2/3*n):n], xlab = 'U', ylab='V', main = '1.3<z<4')



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

