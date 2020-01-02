## ------------------------------------------------------------------------
library(SC19024)
RMonte(FUN=exp,low = 0,up = 5,N=1000,layer = 5, Method = 'control')
RMonte(FUN = function(x){x^2},0,1,1000,5,'control')


## ------------------------------------------------------------------------
# normal monte carlo integration in interval [0,5]
TestVariance(RMonte,100,exp,0,5,1000,1,'None')

# reduce variance by antithetic variable
TestVariance(RMonte,100,exp,0,5,1000,1,'Antithetic')

# reduce variance by control variable
TestVariance(RMonte,100,exp,0,5,1000,1,'control')

#reduce variance by stratification
TestVariance(RMonte,100,exp,0,5,1000,3,'None')

# stratification and control variable at the same time
TestVariance(RMonte,100,exp,0,5,1000,3,'control')


## ------------------------------------------------------------------------
Cboot(DATA = 1:100,sims = 100,Method = 'Perc',alpha = 0.05)

## ----warning=FALSE, message=FALSE----------------------------------------
library(boot)
library(microbenchmark)
microbenchmark(Cboot(DATA = 1:100,sims = 10,Method = 'Perc',alpha = 0.05),boot(1:100,statistic = function(data,index){mean(data[index])},R=10))

## ------------------------------------------------------------------------
library(ggplot2)
NORMAL = C_poly(c(1),c(0,0,-0.5),1,10000,-1000,1000)
EXP = C_poly(c(1),c(0,-1),1,10000,0,1000)
POLY = C_poly(c(0,0,1),c(0),0.1,10000,0,1)


DT = as.data.frame(cbind(NORMAL,EXP,POLY))
ggplot(DT) + geom_histogram(mapping = aes(x=NORMAL,y=..density..),bins = 80) + stat_function(fun=dnorm)
ggplot(DT) + geom_histogram(mapping = aes(x=EXP,y=..density..),bins=80) + stat_function(fun=dexp)
ggplot(DT) + geom_histogram(mapping = aes(x=POLY,y=..density..),bins=80) + stat_function(fun=function(x){3*x^2})

