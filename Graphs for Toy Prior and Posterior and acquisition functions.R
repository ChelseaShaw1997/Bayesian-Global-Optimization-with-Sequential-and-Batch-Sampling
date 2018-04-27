#The work on the prior and posterior distributions follows
#code given in https://www.r-bloggers.com/gaussian-process-regression-with-r/
#by James Keirstead
#This code requires the use of the package 'ggplot2', written by Hadley Wickham and Winston Chang, published in 2016.
#The url for this package is https://cran.r-project.org/web/packages/ggplot2/index.html
#This code requires the use of the package 'plyr', written by Hadley Wickham, published in 2016.
#The url for this package is https://cran.r-project.org/web/packages/plyr/index.html
#This code requires the use of the package 'MASS', written by Brian Ripley et. al., published in 2018.
#The url for this package is https://cran.r-project.org/web/packages/MASS/MASS.pdf
#This code requires the use of the package 'reshape2', written by Hadley Wickham.
#The url for this package is https://cran.r-project.org/web/packages/reshape2/index.html
#This code requires the use of the package 'gridExtra', written by Baptiste Auguie, published in 2017.
#The url for this package is https://cran.r-project.org/web/packages/gridExtra/index.html
#This code requires the use of the package 'mvtnorm', written by Alan Genz et. al., published in 2018.
#The url for this package is https://cran.r-project.org/web/packages/mvtnorm/index.html
require(MASS)
require(mvtnorm)

require(plyr)
require(reshape2)
require(ggplot2)
set.seed(12345)
calcSigma<- function(X1,X2,l=0.4936487) {
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
    }
  }
  return(Sigma)
}

#toy function
quadfunction = function(x){
  ans = -x*(x+1.5)*(x-1.5)*(x-2)
  return(ans)
}


x.star <- seq(-1.8, 2.5, len=50)
sigma <- calcSigma(x.star, x.star)


n.samples <- 3
values <- matrix(rep(0,length(x.star)*n.samples), ncol=n.samples)
for (i in 1:n.samples) {
  # Each column represents a sample from a multivariate normal distribution
  # with zero mean and covariance sigma
  values[,i] <- mvrnorm(1, rep(0, length(x.star)), sigma)
}
values <- cbind(x=x.star,as.data.frame(values))
values <- melt(values,id="x")

#prior graph
fig2a <- ggplot(values,aes(x=x,y=value)) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-2, ymax=2, fill="grey80") +
  geom_line(aes(group=variable, colour = variable)) +
  theme_bw() +
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig2a

#training set
fx = c(-1.25, -0.8, 0, 1.8)
fy = c(quadfunction(-1.15), quadfunction(-0.8), quadfunction(0), quadfunction(1.8))
f = data.frame(x = fx, y = fy)

#posterior calculations
x <- fx
k.xx <- calcSigma(x,x)

k.xxs <- calcSigma(x,x.star)
k.xsx <- calcSigma(x.star,x)
k.xsxs <- calcSigma(x.star,x.star)
f.star.bar <- k.xsx%*%solve(k.xx)%*%fy
cov.f.star <- k.xsxs - k.xsx%*%solve(k.xx)%*%k.xxs

n.samples <- 10000
values <- matrix(rep(0,length(x.star)*n.samples), ncol=n.samples)

for (i in 1:n.samples) {
  values[,i] <- mvrnorm(1, f.star.bar, cov.f.star)
}
values <- cbind(x=x.star,as.data.frame(values))
values <- melt(values,id="x")

#posterior graph
fig2b <- ggplot(values,aes(x=x,y=value)) +
  geom_line(aes(group=variable), colour="grey80") +
  geom_path(data=data.frame(x=x.star,y=f.star.bar), aes(x=x,y=y),colour="red") +
  geom_point(data=f,aes(x=x,y=y)) +
  theme_bw() +
  scale_y_continuous(lim=c(-5,5), name="output, f(x)") +
  xlab("input, x")

fig2b
f.star.bar
x.star

#toy function returning a list
lstquadfunction = function(x){
  ans = -x*(x+1.5)*(x-1.5)*(x-2)
  lst=list(Score = ans, Pred = 0)
  return(lst)
}

variance = diag(cov.f.star)
f.bar.star <- k.xsx%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%f$y
cov.f.star <- k.xsxs - k.xsx%*%solve(k.xx + sigma.n^2*diag(1, ncol(k.xx)))%*%k.xxs

#UCB function
UCB = f.bar.star+2.576*sqrt(variance)
UCBrow1 = seq(-1.8,2.5,length.out = 50)
plottingUCB = cbind(UCBrow1,UCB)
dfplotUCB = as.data.frame(plottingUCB)
colnames(dfplotUCB)=c("x","y")

#plotting the UCB function
UCBplot = ggplot(dfplotUCB) +
  geom_line(aes(x=x,y=y),color="red")
UCBplot
require(gridExtra)

grid.arrange(fig2b,UCBplot)


#function to divide vectors for PI
vectordivision = function(mu,var){
  i = 1
  n = length(mu)
  vec = rep(0,n)
  for(i in 1:n){
    vec[i] = mu[i]/var[i]
    i = i+1
  }
  return(vec)
}


mu.over.sigma = vectordivision(f.star.bar,sqrt(variance))

NormalDistVec = function(input){
  i=1
  n = length(input)
  vec = rep(0,n)
  for(i in 1:n){
    vec[i] = pnorm(input[i])
    i=i+1
  }
  return(vec)
}
densNormVec = function(input){
  i=1
  n = length(input)
  vec = rep(0,n)
  for(i in 1:n){
    vec[i] = dnorm(input[i])
    i=i+1
  }
  return(vec)
}

#PI function
adjmu.over.sigma = vectordivision(f.bar.star-rep(0.01,length(f.bar.star)),sqrt(variance))
PI = NormalDistVec(adjmu.over.sigma)
PIrow1 = seq(-1.8,2.5,length.out=50)
plottingPI = cbind(PIrow1,PI)
dfplotPI = as.data.frame(plottingPI)
colnames(dfplotPI)=c("x","y")
#plotting PI
PIplot = ggplot(dfplotPI) +
  geom_line(aes(x=x,y=y),color="red")
PIplot
grid.arrange(fig2b,PIplot)


#EI
EI = (f.bar.star-rep(1,length(f.bar.star)))*NormalDistVec(adjmu.over.sigma)+sqrt(variance)*densNormVec(adjmu.over.sigma)
EIrow1 = seq(-1.8,2.5,length.out = 50)
plottingEI = cbind(EIrow1,EI)
dfplotEI = as.data.frame(plottingEI)
colnames(dfplotEI)=c("x","y")
#plotting EI
EIplot = ggplot(dfplotEI) +
  geom_line(aes(x=x,y=y),color="red")
EIplot
grid.arrange(fig2b,EIplot)





