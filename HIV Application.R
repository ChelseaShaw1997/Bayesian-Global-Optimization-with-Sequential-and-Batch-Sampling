#The writing of this code was aided by Assistant Proffessor
#Georgios Karagiannis, Durham University
#This code requires use of the package 'deSolve' written by K. Soetaert, T. Petzoldt and R. Woodrow Setzer, published in 2017.
#The url for this package is https://cran.r-project.org/web/packages/deSolve/deSolve.pdf
#This code requires use of the package 'rBayesianOptimization' written by Yachen Yan, published in 2017.
#The url for this package is https://cran.r-project.org/web/packages/rBayesianOptimization/rBayesianOptimization.pdf
#The data values and ODE model were found in "HIV dynamics:modeling, data analysis
#and optimal treatment protocols" by BM Adams et. al., 2005.
library(deSolve)
library(rBayesianOptimization)

#parameters
lam1 = 10000
d1 = 0.01
epsilon = 0
k1=8*(10^(-7))
lam2 = 31.98
d2 = 0.01
f=0.34
k2=10^(-4)
delta = 0.7
m1 = 10^(-5)
m2 = 10^(-5)
NT = 100
c = 13
rho1 = 1
rho2 = 1
lamE = 1
bE = 0.3
Kb = 100
dE = 0.25
Kd = 500
deltaE = 0.1

#the ode model
hiv.ode = function(t,y, parms){
  parms = as.numeric(parms)
  y = as.numeric(y)
  lam1 = parms[20]
  d1 = parms[1]
  epsilon = parms[2]
  k1=parms[18]
  lam2 = parms[21]
  d2 = parms[3]
  f=parms[4]
  k2=parms[19]
  delta = parms[5]
  m1 = parms[6]
  m2 = parms[7]
  NT = parms[8]
  c = parms[9]
  rho1 = parms[10]
  rho2 = parms[11]
  lamE = parms[12]
  bE = parms[13]
  Kb = parms[14]
  dE = parms[15]
  Kd = parms[16]
  deltaE = parms[17]
  dy1 = lam1-(d1*y[1])-(1-epsilon)*k1*y[5]*y[1]
  dy2 = lam2-d2*y[2]-(1-f*epsilon)*k2*y[5]*y[2]
  dy3 = (1-epsilon)*k1*y[5]*y[1]-delta*y[3]-m1*y[6]*y[3]
  dy4 = (1-f*epsilon)*k2*y[5]*y[2]-delta*y[4]-m2*y[6]*y[4]
  dy5 = NT*delta*(y[3]+y[4]) - c*y[5]-((1-epsilon)*rho1*k1*y[1]+(1-f*epsilon)*rho2*k2*y[2])*y[5]
  dy6 = lamE + (bE*(y[3]+y[4])/(y[3]+y[4]+Kb))*y[6]-(dE*(y[3]+y[4])/(y[3]+y[4]+Kd))*y[6]-deltaE*y[6]
  return(list(c(dy1,dy2,dy3,dy4,dy5,dy6)))
}

tvec = seq(t0,tf,1) #time vector
tvec
q=c(0.01,0,0.01,0.34,0.7,10^(-5),10^(-5),100,13,1,1,1,0.3,100,0.25,500,0.1,8*(10^(-7)),10^(-4),10000,31.98)

xbar = ode(y0,tvec,hiv.ode,q)
y40 = xbar[41,-1]
y80 = xbar[81,-1]
y120 = xbar[121,-1]
y160=xbar[161,-1]
y200 = xbar[201,-1]
yexp = as.matrix(t(cbind(y0,y40,y80,y120,y160,y200)))

#this next vector must be changed depending on whih variables are being estimated
fixedq = c(0.01,0,0.01,0.34,0.7,10^(-5),10^(-5),100,13,1,1,1,0.3,100,0.25,500,0.1,8*(10^(-7)),10^(-4),10000)


ode.solver = function(t,varparams){
  params = as.double(c(fixedq,varparams))
  y = ode(y0,t,odename,params)
  return(y)
}
ode.needed.matrix = function(varparams){
  t = tvec
  full.y = ode.solver(t,varparams)
  y0.1 = y0
  y40.1 = full.y[41,-1]
  y80.1 = full.y[81,-1]
  y120.1 = full.y[121,-1]
  y160.1 = full.y[161,-1]
  y200.1 = full.y[201,-1]
  y = as.matrix(t(cbind(y0.1,y40.1,y80.1,y120.1,y160.1,y200.1)))
  return(y)
}
simulated.solution = function(varparams){
  y = as.matrix(ode.needed.matrix(varparams))
  output = matrix(0,6,6)
  i = 1
  w = 1
  for(i in 1:6){
      for(w in 1:6){
        output[i,w]=y[i,w]
        w = w+1
      }
    i = i+1
  }
  return(output)
}

cost.function = function(lam2.var){
  varparams = as.numeric(c(lam2.var))
  t = c(0,40,80,120,160,200)
  simulated.y = simulated.solution(varparams)
  output = -sqrt(mean(simulated.y-yexp)^2)
  return(list("Score" = output, "Pred" = NULL))
}



BGO =BayesianOptimization(FUN = cost.function,
                           bounds = list(lam1.var = c(9500,10500), lam2.var=c(10,50)), init_grid_dt = NULL,
                           init_points = 5, n_iter = 20,
                           acq = "ucb",kernel = list(type = "exponential", power = 2),
                           verbose = TRUE
                           )
BGO2 = BayesianOptimization(FUN = cost.function,
                            bounds = list(k2.var = c(0,0.1),lam1.var = c(9000,11000), lam2.var=c(10,50)), init_grid_dt = NULL,
                            init_points = 5, n_iter = 20,
                            acq = "ucb",kernel = list(type = "exponential", power = 2),
                            verbose = TRUE
)
BGO3 = BayesianOptimization(FUN = cost.function,
                            bounds = list(k1.var = c(0,0.1),k2.var = c(0,0.1),lam1.var = c(9000,11000), lam2.var=c(10,50)), init_grid_dt = NULL,
                            init_points = 5, n_iter = 20,
                            acq = "ucb",kernel = list(type = "exponential", power = 2),
                            verbose = TRUE
)
BGO4 = BayesianOptimization(FUN = cost.function,
                            bounds = list(lam2.var=c(10,50)), init_grid_dt = NULL,
                            init_points = 5, n_iter = 20,
                            acq = "ucb",kernel = list(type = "exponential", power = 2),
                            verbose = TRUE
)


#plotting results

x.1.var = as.vector(BGO4$History[,2])
y.1.var= as.vector(BGO4$History[,3])
plot.1.var = data.frame(x.1.var,y.1.var)
par(mfrow=c(1,1))
plot(plot.1.var, "xlab" = "Value of Lambda 2","ylab"="Value of Cost Function",pch=3)

plot.2.var.3d = matrix(cbind(BGO$History[,2],BGO$History[,3],BGO$History[,4]))
plot.2.var.1 = data.frame(BGO$History[,2],BGO$History[,4])
plot.2.var.2 = data.frame(BGO$History[,3],BGO$History[,4])

par(mfrow = c(1,2))
plot(plot.2.var.1, xlab = "Value of Lambda 1",ylab="Value of Cost Function",pch=3)
plot(plot.2.var.2, xlab = "Value of Lambda 2",ylab="Value of Cost Function",pch=3)

par(mfrow=c(1,3))
plot.3.var.1 = data.frame(BGO2$History[,2],BGO2$History[,5])
plot.3.var.2 = data.frame(BGO2$History[,3],BGO2$History[,5])
plot.3.var.3 = data.frame(BGO2$History[,4],BGO2$History[,5])

plot(plot.3.var.1, xlab = "Value of k2",ylab="Value of Cost Function",pch=3)
plot(plot.3.var.2, xlab = "Value of Lambda 1",ylab="Value of Cost Function",pch=3)
plot(plot.3.var.3, xlab = "Value of Lambda 2",ylab="Value of Cost Function",pch=3)

par(mfrow =c(1,4))
plot.4.var.1 = data.frame(BGO3$History[,2],BGO3$History[,6])
plot.4.var.2 = data.frame(BGO3$History[,3],BGO3$History[,6])
plot.4.var.3 = data.frame(BGO3$History[,4],BGO3$History[,6])
plot.4.var.4 = data.frame(BGO3$History[,5],BGO3$History[,6])

par(mfrow = c(1,4))
plot(plot.4.var.1, xlab = "Value of k1",ylab="Value of Cost Function",pch=3)
plot(plot.4.var.2, xlab = "Value of k2",ylab="Value of Cost Function",pch=3)
plot(plot.4.var.3, xlab = "Value of Lambda 1",ylab="Value of Cost Function",pch=3)
plot(plot.4.var.4, xlab = "Value of Lambda 2",ylab="Value of Cost Function",pch=3)
