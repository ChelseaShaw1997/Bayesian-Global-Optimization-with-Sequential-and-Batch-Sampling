#This code uses the package 'GenSA' written by Gubian, S et. al., published in January 2018. 
#The url for this package is https://cran.r-project.org/web/packages/GenSA/GenSA.pdf

library(GenSA)
library(ggplot2)
#GenSA finds the minimum- here we're maximising, so we need the reverse of the quartic function
minusquarfunction = function(x){
  ans = x*(x+1.5)*(x-1.5)*(x-2)
  return(ans)
}


#Simulated Annealing
SimAn = GenSA(fn = minusquarfunction,lower = -1.8, upper = 2.5, control = list(maxit = 50))

#plotting the graph
dat <- as.data.frame(SimAn$trace.mat)
dat$function.value = -dat$function.value #accounting for maximisation rather than minimisation
dat$current.minimum = -dat$current.minimum #accounting for maximisation rather than minimisation
ggplot(dat, aes(x = nb.steps, y =current.minimum)) +
  geom_line(data = dat, 
            aes(x = nb.steps, y = function.value), 
            position = position_jitter(height = 0, width = 0.05),
            alpha = 0.5)+
  #geom_step(alpha = 0.3) +
  theme_minimal() +
  xlab("Iteration") +
  ylab("Maximum at Current Iteration")
View(dat)
