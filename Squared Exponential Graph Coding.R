#The work on the prior and posterior distributions follows
#code given in https://www.r-bloggers.com/gaussian-process-regression-with-r/
#by James Keirstead

library(GPfit)
require(MASS)


require(plyr)
require(reshape2)
require(ggplot2)
set.seed(12345)
calcSigma0.5 <- function(X1,X2,l=0.5) {
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
    }
  }
  return(Sigma)
}
calcSigma1 <- function(X1,X2,l=1) {
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
    }
  }
  return(Sigma)
}
calcSigma2 <- function(X1,X2,l=2) {
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
    }
  }
  return(Sigma)
}
x.star0.5 <- seq(-1.8, 2.5, len=50)
x.star1 <- seq(-1.8, 2.5, len=50)
x.star2 <- seq(-1.8, 2.5, len=50)
sigma0.5 <- calcSigma0.5(x.star0.5, x.star0.5)
sigma1 <- calcSigma1(x.star1, x.star1)
sigma2 <- calcSigma2(x.star2,x.star2)

fx = c(-1.25, -0.8, 0, 1.8)
fy = c(quadfunction(-1.15), quadfunction(-0.8), quadfunction(0), quadfunction(1.8))
f = data.frame(x = fx, y = fy)

x0.5 <- f$x
x1 <- f$x
x2 <- f$x

k.xx0.5 <- calcSigma0.5(x0.5,x0.5)
k.xxs0.5 <- calcSigma0.5(x0.5,x.star0.5)
k.xsx0.5 <- calcSigma0.5(x.star0.5,x0.5)
k.xsxs0.5 <- calcSigma0.5(x.star0.5,x.star0.5)
f.star.bar0.5 <- k.xsx0.5%*%solve(k.xx0.5)%*%f$y
cov.f.star0.5 <- k.xsxs0.5 - k.xsx0.5%*%solve(k.xx0.5)%*%k.xxs0.5

k.xx1 <- calcSigma1(x1,x1)
k.xxs1 <- calcSigma1(x1,x.star1)
k.xsx1 <- calcSigma1(x.star1,x1)
k.xsxs1 <- calcSigma1(x.star1,x.star1)
f.star.bar1 <- k.xsx1%*%solve(k.xx1)%*%f$y
cov.f.star1 <- k.xsxs1 - k.xsx1%*%solve(k.xx1)%*%k.xxs1

k.xx2<- calcSigma2(x2,x2)
k.xxs2 <- calcSigma2(x2,x.star2)
k.xsx2 <- calcSigma2(x.star2,x2)
k.xsxs2 <- calcSigma2(x.star2,x.star2)
f.star.bar2 <- k.xsx2%*%solve(k.xx2)%*%f$y
cov.f.star2 <- k.xsxs2 - k.xsx2%*%solve(k.xx2)%*%k.xxs2

n.samples0.5 <- 500
n.samples1 <- 500
n.samples2 <- 500
values0.5 <- matrix(rep(0,length(x.star0.5)*n.samples0.5), ncol=n.samples0.5)
values1 <- matrix(rep(0,length(x.star1)*n.samples1), ncol=n.samples1)
values2 <- matrix(rep(0,length(x.star2)*n.samples2), ncol=n.samples2)

i=1
for (i in 1:n.samples0.5) {
  values0.5[,i] <- mvrnorm(1, f.star.bar0.5, cov.f.star0.5)
  i=i+1
}
values0.5 <- cbind(x=x.star0.5,as.data.frame(values0.5))
values0.5 <- melt(values0.5,id="x")

for (i in 1:n.samples1) {
  values1[,i] <- mvrnorm(1, f.star.bar1, cov.f.star1)
  i = i+1
}
values1 <- cbind(x=x.star1,as.data.frame(values1))
values1 <- melt(values1,id="x")

for (i in 1:n.samples2) {
  values2[,i] <- mvrnorm(1, f.star.bar2, cov.f.star2)
  i = i+1
}
values2 <- cbind(x=x.star2,as.data.frame(values2))
values2 <- melt(values2,id="x")

fig2b0.5 <- ggplot(values0.5,aes(x=x,y=value)) +
  geom_line(aes(group=variable), colour="grey80") +
  geom_path(data=data.frame(x=x.star0.5,y=f.star.bar0.5), aes(x=x,y=y),colour="red") +
  geom_point(data=f,aes(x=x,y=y)) +
  theme_bw() +
  scale_y_continuous(lim=c(-5,5), name="output, f(x)") +
  xlab("input, x")

fig2b1 <- ggplot(values1,aes(x=x,y=value)) +
  geom_line(aes(group=variable), colour="grey80") +
  geom_path(data=data.frame(x=x.star1,y=f.star.bar1), aes(x=x,y=y),colour="red") +
  geom_point(data=f,aes(x=x,y=y)) +
  theme_bw() +
  scale_y_continuous(lim=c(-5,5), name="output, f(x)") +
  xlab("input, x")

fig2b2 <- ggplot(values2,aes(x=x,y=value)) +
  geom_line(aes(group=variable), colour="grey80") +
  geom_path(data=data.frame(x=x.star2,y=f.star.bar2), aes(x=x,y=y),colour="red") +
  geom_point(data=f,aes(x=x,y=y)) +
  theme_bw() +
  scale_y_continuous(lim=c(-5,5), name="output, f(x)") +
  xlab("input, x")


require(gridExtra)
grid.arrange(fig2b0.5, fig2b1, fig2b2)


fig2b0.5
fig2b1
fig2b2


