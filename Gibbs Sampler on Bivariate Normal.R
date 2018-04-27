#Guided by code from
#http://www.mas.ncl.ac.uk/~ndjw1/teaching/sim/gibbs/gibbs.html
#which was written by Professor Darren Wilkinson
#This code requires the use of the package 'MASS', written by Brian Ripley et. al., published in 2018.
#The url for this package is https://cran.r-project.org/web/packages/MASS/MASS.pdf
library(MASS)


sig = cbind(c(1,0.95),c(0.95,1))



par(mfrow=c(1,3))
par(mfrow=c(3,2))
bivn = mvrnorm(10000, mu = c(0,0), Sigma = sig)
plot(bivn, col = 1:10000)



X1 = bivn[,1]
X2 = bivn[,2]
hist(X1,50,xlab = "X1")
hist(X2,50, xlab = "X2")

plot(ts(X1), col = 2)
plot(ts(X2), col = 2)

gibbs <- function(x0,y0,n,rho){
  mat = matrix(0,n,2)
  mat[1,]<-c(x0,y0)
  i = 2
  for(i in 2:n){
    w = runif(1,1,2)
    if( w < 1.5){
      x = rnorm(1,rho*mat[i-1,2],sqrt(1-rho^2))
      mat[i,]<-c(x,mat[i-1,2])
    }
    if(w>=1.5){
      y = rnorm(1,rho*mat[i-1,1],sqrt(1-rho^2))
      mat[i,]<-c(mat[i-1,1],y)
    }
    i = i+1
  }
  return(mat)  
}
Gbvn = gibbs(0,0,20000,0.95)
par(mfrow=c(2,2))
plot(Gbvn,col=1:10000, xlab = "Simulated values of X1", ylab = "Simulated values of X2")
plot(Gbvn,type="l")



RSGsimulatedX1 = Gbvn[,1]
RSGsimulatedX2 = Gbvn[,2]
hist(RSGsimulatedX1,50,xlab = "Simulated values of X1")
hist(RSGsimulatedX2,50, xlab = "Simulated Values of X2")

plot(ts(RSGsimulatedX1), col = 4)
plot(ts(RSGsimulatedX2), col = 4)

gibbsseq<-function (n, rho) 
{
  mat <- matrix(ncol = 2, nrow = n)
  x <- 0
  y <- 0
  mat[1, ] <- c(x, y)
  for (i in 2:n) {
    x <- rnorm(1, rho * y, sqrt(1 - rho^2))
    y <- rnorm(1, rho * x, sqrt(1 - rho^2))
    mat[i, ] <- c(x, y)
  }
  mat
}
bvn<-gibbsseq(10000,0.95)
par(mfrow=c(3,2))
plot(bvn,col=1:10000,xlab = "Simulated values of X1", ylab = "Simulated values of X2")
plot(bvn,type="l")
SsimulatedX1 = bvn[,1]
SsimulatedX2 = bvn[,2]
plot(ts(SGsimulatedX1), col = 3)
plot(ts(SGsimulatedX2), col = 3)
SGsimulatedX1 = bvn[,1]
SGsimulatedX2 = bvn[,2]
hist(SGsimulatedX1,50,xlab = "Simulated values of X1")
hist(SGsimulatedX2,50,xlab = "Simulated values of X2")
