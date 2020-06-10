R=6
mu=c(-9, -5.4, -1.8, 1.8, 5.4, 9)
sigmar=rep(1/sqrt(10), R)
pi=rep(1/R, R)
alpha=3/4
beta=1/2
theta=20

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

#test 
X = c(rnorm(1000,mean=-9,sd=1/sqrt(10)),
      rnorm(1000,mean=-5.4,sd=1/sqrt(10)),
      rnorm(1000,mean=-1.8,sd=1/sqrt(10)),
      rnorm(1000,mean=1.8,sd=1/sqrt(10)),
      rnorm(1000,mean=5.4,sd=1/sqrt(10)),
      rnorm(1000,mean=9,sd=1/sqrt(10)))

#Sample quantile
quantile(X,.95)

#Quantile from F_inv
F_inv(.95,pi ,mu, sigmar)

V=runif(1000)
W=runif(1000)


C <- function(u) {
  (1-beta)*u^(1-alpha)*v^(-beta)*(u^(-theta*alpha)+v^(-theta*beta)-1)^(-1/theta)+
    beta*u^(1-alpha)*v^(-beta*(1+theta))*(u^(-theta*alpha)+v^(-theta*beta)-1)^(-1/theta-1)-w
}
U=rep(0,1000)
for (i in 1:1000) {
  v=V[i]
  w=W[i]
  U[i]=uniroot(C, lower = 0, upper = 1)$root
}

X1=sapply(V, F_inv, w=pi, u=mu, s=sigmar)
X2=sapply(U, F_inv, w=pi, u=mu, s=sigmar)


copula=data.frame(V, U)
# write.csv(data, "UV.csv", row.names = FALSE)

# Z=qnorm(U);
# W=qnorm(V);
# data=data.frame(Z, W)
# plot(Z,W)
# write.csv(data, "ZW.csv", row.names = FALSE)

library("mixtools")
est=mvnormalmixEM(data, lambda = NULL, mu = NULL, sigma = NULL, k = 3,
              arbmean = TRUE, arbvar = TRUE, epsilon = 1e-08, 
              maxit = 10000, verb = FALSE)

est$lambda
est$mu
est$sigma

x.1<-rmvnorm(round(est$lambda[1],2)*1000, est$mu[[1]], est$sigma[[1]])
x.2 <- rmvnorm(round(est$lambda[2],2)*1000, est$mu[[2]], est$sigma[[2]])
x.3 <-rmvnorm(round(est$lambda[3],2)*1000, est$mu[[3]], est$sigma[[3]])

X.1 <- rbind(x.1, x.2, x.3)
pred=data.frame(X.1)
colnames(pred)=c("Z", "W")
copula.pred=data.frame(pnorm(pred$W), pnorm(pred$Z))
colnames(copula.pred)=c("V", "U")

ggplot(copula) +
  geom_point(aes(x=V,y=U, color='blue'))+
  geom_point(data=copula.pred,aes(x=V,y=U, color = 'red')) +
  scale_color_discrete(name='', labels=c('Pred','Raw'))

par(mfrow=c(1,1))
plot(data)
par(new=TRUE)
points(pred, col='green')

library("ggplot2")
ggplot(data) +
  geom_point(aes(x=Z,y=W, color='blue'))+
  geom_point(data=pred,aes(x=Z,y=W, color = 'red')) +
  scale_color_discrete(name='', labels=c('Pred','Raw'))

