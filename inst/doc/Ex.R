## ----warning=FALSE, message=FALSE----------------------------------------
x = rnorm(1000) #simulate 100 data
x_mean = mean(x)
x_var = var(x)
text1 = paste('mean:',x_mean, sep = '')
text2 = paste('var:', x_var, sep = '')
print(text1)
print(text2)

## ----warning=FALSE, message=FALSE----------------------------------------
f <- function(x){
  return(2^(-3/2)*pi^(-0.25)*x^0.5*exp(-x/2))
}
x = seq(from=0, to=5, by=0.01) #cut the interval
y = f(x)
plot(x,y,type='l') #plot

## ----warning=FALSE, message=FALSE----------------------------------------
data('cars') #load dataset
library(knitr)
out = lm(dist~speed,cars)#fit linear model
kable(summary(out)$coefficients)

## ----warning=F, message=F------------------------------------------------
#genereate U(0,1), set seed = 357
set.seed(357)
U = runif(1000)
#inverse tranformation
U_sp = vector(mode='numeric', length=1000)
U_sp[U < 0.1] = 0
U_sp[U < 0.3 & U >= 0.1] = 1
U_sp[U < 0.5 & U >= 0.3] = 2
U_sp[U < 0.7 & U >= 0.5] = 3
U_sp[U >= 0.7] = 4
freq <- table(U_sp)
# add a row of theatrical result
gene_prob = freq/1000
theory_prob = c(0.1,0.2,0.2,0.2,0.3)
tb = rbind(freq,gene_prob, theory_prob)
knitr::kable(tb)
#generate with Rsample
U_real = sample(c(0,1,2,3,4),replace = T, size = 1000, prob = c(0.1,0.2,0.2,0.2,0.3))
table(U_real) -> freq2
gene_prob2 = freq2/1000
tb2 = rbind(freq2, gene_prob2, theory_prob)
knitr::kable(tb2)

## ----warning=F, message=F------------------------------------------------
# write a function
BETA <- function(a,b,n){

  #let c = beta(a,b)
  #compute c
  result = vector()
  #generate random variables until we have enough samples(n)
  while(TRUE){
    c = factorial(a+b-1) / factorial(a-1) / factorial(b - 1)
  #generate variables
    len = n*c*2 
    U = runif(len)
    G0 = runif(len)
    F0 = c*G0^(a-1)*(1-G0)^(b-1)
    result_temp = G0[(F0) >= c*U]
    result = c(result, result_temp)
    if(length(result)>=n){
      break;
    }
  }
  return(result[1:n])
}

result = BETA(3,2,1000)
#draw the density curve 
x = seq(0,1,0.01)
y = (x)^2*(1-x)^1/beta(3,2)
par(new=T)
hist(result, breaks=30,xlab = '', probability = T)
curve(dbeta(x,3,2),add=TRUE, lwd=2)

## ----warning=F, message=F------------------------------------------------
n = 1000 # number of sample
r = 4 
beta = 2 

Lambda = rgamma(n,r,beta)

Y = rexp(n,Lambda)

hist(Y)

## ----warning=FALSE,message=FALSE-----------------------------------------
library(knitr)
library(ggplot2)
UP = pi/3
U = runif(10000, 0, UP)
G = sin(U)
CONTROL = U - U^3/6
ECONTROL = pi/6 - pi^3/648
L = lm(G ~ CONTROL)
cstar = L$coefficients[2]
G_adjust = G - cstar*(CONTROL - ECONTROL)
result_adj = mean(G_adjust)*UP
result_0 = mean(G)*UP
var_0 = var(G)*UP^2
var_adjust = var(G_adjust)*UP^2
percentage_empirical = 1 - var_adjust/var_0
percentage_theorical = summary(L)$r.squared
result = as.matrix(c(result_0,result_adj))
result = cbind(result,c(sqrt(var_0),sqrt(var_adjust)))
result = cbind(result, c(0.5,0.5))
rownames(result) <- c('unadjust','adjust')
colnames(result) = c('estimate','standard error', 'True')
kable(result)
100*c(percentage_empirical,percentage_theorical) ->tb
names(tb) <-c('empirical', 'theorical')
tb = as.matrix(tb)
colnames(tb) <- c('percentage of variance reduced(control variate)')
kable(tb)

## ----warning=FALSE, message=FALSE----------------------------------------
# up is the upper bound(pi/3)
x = seq(0,UP,0.01)
y = sin(x)
data = as.data.frame(cbind(x,y))
g = ggplot(data, aes(x,y))
g + geom_line() + scale_x_continuous(breaks = seq(0,1.5,0.1)) + scale_y_continuous(limits = c(0,1),breaks = seq(0,1.5,0.1))

## ----warning=FALSE,message=FALSE-----------------------------------------
y2 = 18/pi^2*x
y3 = y/y2
data = as.data.frame(cbind(x,y,y2,y3))
g = ggplot(data)
g + geom_line(aes(x,y)) + geom_line(aes(x,y2), color='green')+ geom_line(aes(x,y3), color='red')

## ----warning=FALSE, message=FALSE----------------------------------------
rf <- function(n){
  C = pi^2/9
  U = runif(n)
  return(sqrt(U*C))
}
X = rf(10000)
G_imp = sin(X)
f_imp = 18/pi^2*X
result_imp = mean(G_imp/f_imp)
var_imp = var(G_imp/f_imp)
result_imp = as.matrix(c(result_0,result_imp))
result_imp = cbind(result_imp,c(sqrt(var_0),sqrt(var_imp)))
result_imp = cbind(result_imp, c(0.5,0.5))
rownames(result_imp) <- c('unadjust','adjust')
colnames(result_imp) = c('estimate','standard error', 'True')
kable(result_imp)
rm(list = ls())#rmove all the data and variables in this question

## ----warning=FALSE, message=FALSE----------------------------------------
X = runif(10000)
Y = 1-X
G1 = exp(-X)/(1+X^2)
G2 = exp(-Y)/(1+Y^2)
G = (G1+G2)/2
X2 = runif(10000)
G3 = exp(-X2)/(1+X2^2)
result_original = c(mean(c(G1,G3)), sd(c(G1,G3)))
result_anti = c(mean(G), sd(G))    
result =matrix(c(result_anti,result_original),2,2,byrow = T)
colnames(result) <- c('estimate', 'standard error')
rownames(result) <- c('anti', 'normal')
percentage_reduced = 1-var(G)/var(c(G1,G3))
kable(result)
print(percentage_reduced)

## ----warning=FALSE, message=FALSE----------------------------------------
#function of cdf X
pdf <- function(x){
  return(exp(-x)/(1-exp(-1)))
}
cdf <- function(x){
  return((1-exp(-x))/(1-exp(-1)))
}
# function to generate random variable X folloing distribution f
rf <- function(n){
  U <- runif(n)
  return(-log(1-(1-exp(-1))*U))
}
# original method
X = rf(10000)
G = exp(-X)/(1+X^2)
f = pdf(X)
result = mean(G/f)
v = var(G/f)
# stratified importance sampling
X0 = rf(20000)
f_s = G_s = list(5)
X_s = list(5)
result_s = var_s = vector(mode = 'numeric',length = 5)
X_s[[1]] = X0[cdf(X0)<=0.2][1:2000]
X_s[[2]] = X0[cdf(X0)>0.2&cdf(X0)<=0.4][1:2000]
X_s[[3]] = X0[cdf(X0)>0.4&cdf(X0)<=0.6][1:2000]
X_s[[4]] = X0[cdf(X0)>0.6&cdf(X0)<=0.8][1:2000]
X_s[[5]] = X0[cdf(X0)>0.8&cdf(X0)<=1][1:2000]
for(i in 1:5){
  G_s[[i]] = exp(-X_s[[i]])/(1+X_s[[i]]^2)
  f_s[[i]] = pdf(X_s[[i]])*5
  result_s[i] = mean(G_s[[i]]/f_s[[i]])
  var_s[i] = var(G_s[[i]]/f_s[[i]])
}
est_s  = sum(result_s)
v_s = sum(var_s)
est_importance = c(result,sqrt(v))
est_stratify_importance = c(est_s,sqrt(v_s))
mt = as.matrix(rbind(est_importance,est_stratify_importance))
colnames(mt) = c('est', 'sd')
kable(mt)

## ----warning=FALSE, message=FALSE----------------------------------------
set.seed(4257)
library(ggplot2)
# get 1000 times of Chisq2(2) Sample and compute the confidence interval
generate = function(n, alpha){
  X = rchisq(n,2)
  B = abs(sqrt(var(X)/n)*qt(1-alpha/2, df=n-1))
  return(c(mean(X)-B,mean(X)+B))
}
RESULT = replicate(1000, expr = generate(20,0.05))
REALmean = 2
temp1 = RESULT[1,]<2
temp2 = RESULT[2,]>2
Proportion = sum(temp1&temp2)/1000 * 100# the proportion that confidence interval contains the real mean
knitr::kable(Proportion)
#this is the graph

## ----warning=FALSE, message=FALSE----------------------------------------
data_0 = as.data.frame(cbind(t(RESULT),1:1000))
data_0$fac = as.factor(temp1&temp2)
ggplot(data_0) + geom_linerange(aes(x=V3, ymin=V1,ymax = V2,color=fac),size = 1) + scale_color_manual(values=c("#999999", "#E69F00"))

## ----warning=FALSE, message=FALSE----------------------------------------
getskewness = function(N){
  #N is the sample size of the N(0,1) distribution
  X = rnorm(N)
  m = mean(X)
  A = sum((X - m)^3)/N
  B = (sum((X-m)^2)/N)^1.5
  return (A/B)
}
n = 2000
m = 10000
#generate the samples
Y = replicate(m, expr = getskewness(n))
# f(x) in F(x) = p, using normal distribution as an approximation
sigma = 6*(n-2)/(n+1)/(n+3)
p = c(0.025,0.05,0.95,0.975)
Z = sort(Y)
place_est = Z[p*m]
density_est = dnorm(place_est, 0, sqrt(sigma))
place_sd = sqrt(p*(1-p)/m/density_est^2)
place_limit = qnorm(p,0,sqrt(6/n))
#draw the table
result = rbind(place_est,place_limit,place_sd)
colnames(result) = c('0.025','0.05','0.95','0.975')
rownames(result) = c('Monte Carlo estimation', 'Limit estimation','se for Monte Carlo estimation' )
knitr::kable(result)

## ----warning=FALSE, message=FALSE----------------------------------------
data0 = as.data.frame(Y)
X = seq(-0.5,0.5,0.001)
N_X = dnorm(X,0,sqrt(6/n)) 
data1 = as.data.frame(cbind(X,N_X))
ggplot() + geom_histogram(data = data0,mapping = aes(x = Y,y=..density..),bins = 100) + geom_line(data = data1, aes(X,N_X, color = 'red'))

## ----warning=FALSE, message=FALSE----------------------------------------
library(ggplot2)
#cunction to compute skewness
  skew <- function(X){
    M = mean(X)
    m3 <- mean((X - M)^3)
    m2 <- mean((X - M)^2)
    return(m3/m2^1.5)
  }
#beta distribution
POW <- function(alpha){
  m = 10000
  n = 100
  threshold = 0.05
  DATA = list()
  SK = numeric(m)
  cv = qnorm(1-threshold/2, 0, sqrt(6*(n-2)/(n+1)/(n+3)))
  for(i in 1:m){
    DATA[[i]] = rbeta(n,alpha,alpha)
  }
  SK = sapply(DATA, skew)
  return(list(sum(abs(SK)>=cv)/m,SK))
}
#t distribution
POW2 <- function(miu){
  m = 10000
  n = 100
  threshold = 0.05
  DATA = list()
  SK = numeric(m)
  cv = qnorm(1-threshold/2, 0, sqrt(6*(n-2)/(n+1)/(n+3)))
  for(i in 1:m){
    DATA[[i]] = rt(n,miu)
  }
  SK = sapply(DATA, skew)
  return(list(sum(abs(SK)>=cv)/m,SK))
}
#normal distribution
POW_N <- function(){
  m = 10000
  n = 100
  threshold = 0.05
  DATA = list()
  SK = numeric(m)
  cv = qnorm(1-threshold/2, 0, sqrt(6*(n-2)/(n+1)/(n+3)))
  for(i in 1:m){
    DATA[[i]] = rnorm(n)
  }
  SK = sapply(DATA, skew)
  return(list(sum(abs(SK)>=cv)/m,SK))
}

## ----warning=FALSE, message=FALSE----------------------------------------
p = numeric(6)
p[1] = POW(1)[[1]]
p[2] = POW(10)[[1]]
p[3] = POW(20)[[1]]
p[4] = POW2(2)[[1]]
p[5] = POW2(5)[[1]]
p[6] = POW2(10)[[1]]
names(p) = c('beta1','beta10','beta20','t2','t5','t10')
knitr::kable(p)
#plot the density of normal,beta,t distribution
x = seq(-5,5,0.05)
y_norm = dnorm(x)
y_1 = dbeta(x,1,1)
y_2 = dbeta(x,10,10)
y_3 = dbeta(x,20,20)
y_4 = dt(x,2)
y_5 = dt(x,5)
y_6 = dt(x,10)
DT = as.data.frame(cbind(x,y_norm,y_1,y_2,y_3,y_4,y_5,y_6))
ggplot(data=DT) + geom_line(mapping = aes(x,y_1))+ geom_line(mapping = aes(x,y_4)) + geom_line(aes(x,y_norm),color='red')
SK_N = POW_N()[[2]]
SK_beta_1 = POW(1)[[2]]
SK_t_2 = POW2(2)[[2]]
SK = as.data.frame(cbind(SK_N,SK_t_2,SK_beta_1))

ggplot(data = SK) + geom_density(mapping = aes(x=SK_N),color='red')+ geom_density(mapping = aes(x=SK_beta_1),color='blue')+ geom_density(mapping = aes(x=SK_t_2))

## ----warning=FALSE, message=FALSE----------------------------------------
#p is the probility of the beta distribution or t distribution
POW2_1 <- function(alpha,p){
  m = 10000
  n = 100
  threshold = 0.05
  DATA = list()
  SK = numeric(m)
  cv = qnorm(1-threshold/2, 0, sqrt(6*(n-2)/(n+1)/(n+3)))
  for(i in 1:m){
    flag = sample(c(0,1),replace=TRUE,size=n,prob = c(1-p,p))
    DATA[[i]] = rbeta(n,alpha,alpha)*flag + rnorm(n)*(1-flag)
  }
  SK = sapply(DATA, skew)
  return(list(sum(abs(SK)>=cv)/m,SK))
}
POW2_2 <- function(miu,p){
  m = 10000
  n = 100
  threshold = 0.05
  DATA = list()
  SK = numeric(m)
  cv = qnorm(1-threshold/2, 0, sqrt(6*(n-2)/(n+1)/(n+3)))
  for(i in 1:m){
    flag = sample(c(0,1),replace=TRUE,size=n,prob = c(1-p,p))
    DATA[[i]] = rt(n,miu)*flag + rnorm(n)*(1-flag)
  }
  SK = sapply(DATA, skew)
  return(list(sum(abs(SK)>=cv)/m,SK))
}
#p=0.01,0.1,0.5
#alpha = 1,10
#miu = 2,5
result = list()
result[[1]] = c(POW2_1(1,0.01)[[1]],POW2_1(10,0.01)[[1]],POW2_2(2,0.01)[[1]],POW2_2(5,0.01)[[1]])
result[[2]] = c(POW2_1(1,0.1)[[1]],POW2_1(10,0.1)[[1]],POW2_2(2,0.1)[[1]],POW2_2(5,0.1)[[1]])
result[[3]] = c(POW2_1(1,0.5)[[1]],POW2_1(10,0.5)[[1]],POW2_2(2,0.5)[[1]],POW2_2(5,0.5)[[1]])
names(result) = c('proportion-0.01','proportion-0.1','proportion-0.5')
r = as.data.frame(result)
rownames(r) = c('alpha1','alpha10','miu2','miu5')
knitr::kable(r)

## ----warning=FALSE, message=FALSE----------------------------------------
n = 20
Data1 = numeric(10000)
Data2 = numeric(10000)
Data3 = numeric(10000)
Data4 = numeric(10000)
for(i in 1:10000){
Data1[i] = t.test(rchisq(n,1),mu=1,alternative = 'two.sided')$p.value
Data2[i] = t.test(runif(n,0,2),mu=1,alternative = 'two.sided')$p.value
Data3[i] = t.test(rexp(n,1),mu=1,alternative = 'two.sided')$p.value
Data4[i] = t.test(rnorm(n,1,1),mu=1,alternative = 'two.sided')$p.value
}
p1 = sum(abs(Data1)<0.05)/10000
p2 = sum(abs(Data2)<0.05)/10000
p3 = sum(abs(Data3)<0.05)/10000
p4 = sum(abs(Data4)<0.05)/10000
result = c(p1,p2,p3,p4)
names(result) = c('Chi1','U02','exp1','norm1')
knitr::kable(result)
#y1 = t.test(X1, )

## ----warning=FALSE, message=FALSE----------------------------------------
library(ggplot2)
library(bootstrap)
library(corrplot)
data(scor)
ggplot(data=scor,mapping = aes(x=mec,y=vec))+geom_point()
ggplot(data=scor,mapping = aes(x=mec,y=alg))+geom_point()
ggplot(data=scor,mapping = aes(x=mec,y=alg))+geom_point()
ggplot(data=scor,mapping = aes(x=mec,y=ana))+geom_point()
ggplot(data=scor,mapping = aes(x=mec,y=sta))+geom_point()
ggplot(data=scor,mapping = aes(x=vec,y=alg))+geom_point()
ggplot(data=scor,mapping = aes(x=vec,y=ana))+geom_point()
ggplot(data=scor,mapping = aes(x=vec,y=sta))+geom_point()
ggplot(data=scor,mapping = aes(x=alg,y=ana))+geom_point()
ggplot(data=scor,mapping = aes(x=alg,y=sta))+geom_point()
ggplot(data=scor,mapping = aes(x=ana,y=sta))+geom_point()
knitr::kable(cor(scor))
corrplot(cor(scor))

## ----warning=FALSE, message=FALSE----------------------------------------
#bootstrap
library(boot)
r_12 = function(x,i){
  return(cor(x$mec[i],x$vec[i]))
}
r_34 = function(x,i){
  return(cor(x$alg[i],x$ana[i]))
}
r_35 = function(x,i){
  return(cor(x$alg[i],x$sta[i]))
}
r_45 = function(x,i){
  return(cor(x$ana[i],x$sta[i]))
}
result_12 = boot(scor,r_12,1000)
result_34 = boot(scor,r_34,1000)
result_35 = boot(scor,r_35,1000)
result_45 = boot(scor,r_45,1000)
result = c(sd(result_12$t),sd(result_34$t),sd(result_35$t),sd(result_45$t))
names(result) = c('12','34','35','45')
knitr::kable(result)

## ----warning=FALSE, message=FALSE----------------------------------------
data_norm = matrix(rnorm(50*500),ncol=50,nrow=500)
data_chisq = matrix(rchisq(50*500,5),ncol=50,nrow=500)
interval_norm = list()
interval_chisq = list()
SK = function(x,i){
  EX = mean(x[i])
  EX2 = sum((x[i] - EX)^2)/length(i)
  EX3 = sum((x[i] - EX)^3)/length(i)
  return(EX3/EX2^1.5)
}
#do the boostrap to estimate sd.
for(i in 1:500){
  bnorm = boot(data_norm[i,],SK,500)
  bchisq = boot(data_chisq[i,],SK,500)
  interval_norm[[i]] = c(boot.ci(bnorm,conf = 0.95,type = 'norm'),boot.ci(bnorm,conf = 0.95,type = 'basic'),boot.ci(bnorm,conf = 0.95,type = 'perc'))
    interval_chisq[[i]] = c(boot.ci(bchisq,conf = 0.95,type = 'norm'),boot.ci(bchisq,conf = 0.95,type = 'basic'),boot.ci(bchisq,conf = 0.95,type = 'perc'))
}
real_norm = 0
real_chisq = sqrt(8/5)
# two functions that compute the samples in the intervals and fall in the left/right of the interval
get_in_norm = function(x){
  in_normal = (x$normal[2]<0)&(x$normal[3]>0)
  left_normal = x$normal[2]>0
  right_normal = x$normal[3]<0
  in_basic = (x$basic[4]<0)&(x$basic[5]>0)
  left_basic = x$basic[4]>0
  right_basic = x$basic[5]<0
  in_perc = (x$perc[4]<0)&(x$perc[5]>0)
  left_perc = x$perc[4]>0
  right_perc = x$perc[5]<0
  return(c(in_normal,left_normal, right_normal,in_basic,left_basic,right_basic,in_perc,left_perc,right_perc))
}
get_in_chisq = function(x){
  in_normal = (x$normal[2]<real_chisq)&(x$normal[3]>real_chisq)
  left_normal = x$normal[2]>real_chisq
  right_normal = x$normal[3]<real_chisq
  in_basic = (x$basic[4]<real_chisq)&(x$basic[5]>real_chisq)
  left_basic = x$basic[4]>real_chisq
  right_basic = x$basic[5]<real_chisq
  in_perc = (x$perc[4]<real_chisq)&(x$perc[5]>real_chisq)
  left_perc = x$perc[4]>real_chisq
  right_perc = x$perc[5]<real_chisq
  return(c(in_normal,left_normal, right_normal,in_basic,left_basic,right_basic,in_perc,left_perc,right_perc))
}
percent_normal = sapply(interval_norm, get_in_norm)
percent_chisq = sapply(interval_chisq, get_in_chisq)
result = matrix(ncol = 9,nrow = 2)
for(i in 1:9){
  result[1,i] = mean(percent_normal[i,])
  result[2,i] = mean(percent_chisq[i,])
}
colnames(result) = c('in_normal','left_normal','right_normal',
   'in_basic','left_basic','right_basic',
   'in_perc','left_perc','right_perc')
rownames(result) = c('norm','chisq5')
knitr::kable(result)

## ----warning=FALSE, message=FALSE----------------------------------------
library(bootstrap)
library(boot)
data(scor)
Length = nrow(scor)
#compute the estiate theta-hat.since cov function in R have been adjusted, I use (n-1)/n*cov as the MLE of covariance matrix estimation
Sigma_hat = cov(scor)*(Length-1)/Length
EIGEN = eigen(Sigma_hat)
theta_hat = EIGEN$values[1]/sum(EIGEN$values)
#do the jacknife 
JACK = vector(length=Length)
for(i in 1:Length){
  Sigma_hat = cov(scor[-i,])*(Length-2)/(Length-1)
  EIGEN = eigen(Sigma_hat)
  JACK[i] = EIGEN$values[1]/sum(EIGEN$values)
}
BIAS = (mean(JACK)-theta_hat)*(Length-1)
SE = sqrt((Length-1)/(Length)*sum((JACK-mean(JACK))^2))
##########################################################
#Another solution using jacknife function in boot package 
SCOR = as.matrix(scor)
est = function(X){
  Sigma_hat = cov(SCOR[X,])*(Length-2)/(Length-1)
  EIGEN = eigen(Sigma_hat)
  return = EIGEN$values[1]/sum(EIGEN$values)
}
RESULT = jackknife(1:Length,est)

##########################################################

result = c(theta_hat,BIAS,SE)
names(result) = c('theta_hat','BIAS','SE')
knitr::kable(result)

## ----warning=FALSE, message=FALSE----------------------------------------
library(DAAG)
attach(ironslag)
L10 <- lm(magnetic ~ chemical)
L20 <- lm(magnetic ~ chemical + I(chemical^2))
L30 <- lm(log(magnetic) ~ chemical)
L40 <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3))
Length = nrow(ironslag)
result = matrix(ncol = 4,nrow = Length)
for(i in 1:Length){
  mag = magnetic[-i]
  chem = chemical[-i]
  L1 <- lm(mag ~ chem)
  L2 <- lm(mag ~ chem + I(chem^2))
  L3 <- lm(log(mag) ~ chem)
  L4 <- lm(mag ~ chem + I(chem^2) + I(chem^3))
  result[i,1] = magnetic[i] - predict(L1,data.frame(chem=chemical[i])) 
  result[i,2] = magnetic[i] - predict(L2,data.frame(chem=chemical[i])) 
  result[i,3] = magnetic[i] - exp(predict(L3,data.frame(chem=chemical[i])))
  result[i,4] = magnetic[i] - predict(L4,data.frame(chem=chemical[i])) 
}
CR = c(mean(result[,1]^2),mean(result[,2]^2),mean(result[,3]^2),mean(result[,4]^2))
R_square = c(summary(L10)$adj.r.squared,summary(L20)$adj.r.squared,summary(L30)$adj.r.squared,summary(L40)$adj.r.squared)
names(CR) = c('LM','quadratic','exp','cubic')
names(R_square) = c('LM','quadratic','exp','cubic')

## ------------------------------------------------------------------------
knitr::kable(CR)

## ------------------------------------------------------------------------
knitr::kable(R_square)

## ------------------------------------------------------------------------
library(ggplot2)
MYdata = ironslag
MYdata$L1 = L10$fitted.values
MYdata$L2 = L20$fitted.values
MYdata$L3 = L30$fitted.values
MYdata$L4 = L40$fitted.values
ggplot(data = MYdata)+geom_point(mapping = aes(x=chemical,y=magnetic)) + geom_line(mapping = aes(x=chemical,y=L1),color='red')+geom_line(mapping = aes(x=chemical,y=L2),color='yellow')+geom_line(mapping = aes(x=chemical,y=exp(L3)),color='green')+geom_line(mapping = aes(x=chemical,y=L4),color='blue')

## ----warning=FALSE,message=FALSE,cache=FALSE-----------------------------
library(ggplot2)
n1 = 20
n2 = 30
m = 1000
result = replicate(m, expr={
  X = rnorm(n1)
  Y = rnorm(n2)
  L = c(X,Y)
  x = X-mean(X)
  y = Y-mean(Y)
  L = c(x,y)
  count5 = max(sum(x > max(y)) + sum(x < min(y)),sum(y > max(x)) + sum(y < min(x)))
  temp = numeric(length = 999)
  for (i in 1:999){
    permu = sample(1:(n1+n2),n1+n2)
    X2 = L[permu[1:n1]]
    Y2 = L[permu[(n1+1):(n1+n2)]]
    x = X2 - mean(X2)
    y = Y2 - mean(Y2)
    temp[i] = max(sum(x > max(y)) + sum(x < min(y)),sum(y > max(x)) + sum(y < min(x)))
  }
  sum(count5>temp)
})
pt = as.data.frame(result/1000)
ggplot(data = pt)+geom_histogram(mapping = aes(result/1000,..density..),bins = 20)
print(sum(pt<=0.95))

## ----warning=FALSE,message=FALSE-----------------------------------------
library(MASS)
library(Ball)
library(energy)
library(boot)
library(ggplot2)
N = 100
getpower = function(N){
POWER = matrix(0,ncol = 4,nrow = 100)
  for(i in 1:100){
X = mvrnorm(N,c(0,0),matrix(c(1,0,0,1),nrow=2,byrow = T))
E = mvrnorm(N,c(0,0),matrix(c(1,0,0,1),nrow=2,byrow = T))
Y1 = X/4+E
Y2 = X/4*E
POWER[i,1] = bcov.test(X,Y1,R=999)$p.value
POWER[i,2] = bcov.test(X,Y2,R=999)$p.value
POWER[i,3] = dcor.test(X,Y1,R=999)$p.value
POWER[i,4] = dcor.test(X,Y2,R=999)$p.value
  }
return (colSums(POWER<0.05)/100)
}
POWER = matrix(0,ncol=4,nrow=80)
for(i in 1:40){
  POWER[i,] = getpower(i)
}
DT = as.data.frame(cbind(1:80,POWER))
colnames(DT) = c('n','B1','B2','D1','D2')
#ball, Y1-X
ggplot(data=DT) + geom_point(mapping = aes(n,B1))
#distance, Y1-X
ggplot(data=DT) + geom_point(mapping = aes(n,D1))
#ball,Y2-X
ggplot(data=DT) + geom_point(mapping = aes(n,B2))
#distance,Y2-x
ggplot(data=DT) + geom_point(mapping = aes(n,D2))

## ----warning=FALSE,message=FALSE-----------------------------------------
library(ggplot2)
sigma = 1
L = 10000
Series = numeric(L)
Uniform = runif(L)
k=0
for (i in 2:L){
  y = rnorm(1,Series[i-1],sigma)
  if(Uniform[i]<=exp(-abs(y)+abs(Series[i-1]))){
    Series[i] = y
    k=k+1
  }
  else{
    Series[i] = Series[i-1]
  }
}
AC_rate = k/L
print(AC_rate)
DATA_0 = as.data.frame(cbind(Series,1:L))
LP = function(x){1/2*exp(-abs(x))}
ggplot(data = DATA_0)+geom_histogram(mapping = aes(Series,..density..),bins = 70,fill='#FFCC00')+stat_function(fun=LP)+labs(title = "sigma=1")+theme(plot.title=element_text(hjust=0.5))+annotate(geom = 'text',x=5,y=0.4,label=as.character(AC_rate))
DATA_last = as.data.frame(cbind(Series[9000:10000],9000:10000))
ggplot(data = DATA_0)+geom_line(aes(V2,Series))+labs(title = "all steps")
ggplot(data=DATA_last)+geom_line(aes(V2,V1))+labs(title = "last 1000 steps")

## ----warning=FALSE,message=FALSE-----------------------------------------
library(ggplot2)
sigma = 2
L = 10000
Series = numeric(L)
Uniform = runif(L)
k=0
for (i in 2:L){
  y = rnorm(1,Series[i-1],sigma)
  if(Uniform[i]<=exp(-abs(y)+abs(Series[i-1]))){
    Series[i] = y
    k=k+1
  }
  else{
    Series[i] = Series[i-1]
  }
}
AC_rate = k/L
print(AC_rate)
DATA_0 = as.data.frame(cbind(Series,1:L))
LP = function(x){1/2*exp(-abs(x))}
ggplot(data = DATA_0)+geom_histogram(mapping = aes(Series,..density..),bins = 70,fill='#FFCC00')+stat_function(fun=LP)+labs(title = "sigma=2")+theme(plot.title=element_text(hjust=0.5))+annotate(geom = 'text',x=5,y=0.4,label=as.character(AC_rate))
DATA_last = as.data.frame(cbind(Series[9000:10000],9000:10000))
ggplot(data = DATA_0)+geom_line(aes(V2,Series))+labs(title = "all steps")
ggplot(data=DATA_last)+geom_line(aes(V2,V1))+labs(title = "last 1000 steps")

## ----warning=FALSE,message=FALSE-----------------------------------------
library(ggplot2)
sigma = 5
L = 10000
Series = numeric(L)
Uniform = runif(L)
k=0
for (i in 2:L){
  y = rnorm(1,Series[i-1],sigma)
  if(Uniform[i]<=exp(-abs(y)+abs(Series[i-1]))){
    Series[i] = y
    k=k+1
  }
  else{
    Series[i] = Series[i-1]
  }
}
AC_rate = k/L
print(AC_rate)
DATA_0 = as.data.frame(cbind(Series,1:L))
LP = function(x){1/2*exp(-abs(x))}
ggplot(data = DATA_0)+geom_histogram(mapping = aes(Series,..density..),bins = 70,fill='#FFCC00')+stat_function(fun=LP)+labs(title = "sigma=5")+theme(plot.title=element_text(hjust=0.5))+annotate(geom = 'text',x=5,y=0.4,label=as.character(AC_rate))
DATA_last = as.data.frame(cbind(Series[9000:10000],9000:10000))
ggplot(data = DATA_0)+geom_line(aes(V2,Series))+labs(title = "all steps")
ggplot(data=DATA_last)+geom_line(aes(V2,V1))+labs(title = "last 1000 steps")

## ----warning=FALSE,message=FALSE-----------------------------------------
library(ggplot2)
sigma = 10
L = 10000
Series = numeric(L)
Uniform = runif(L)
k=0
for (i in 2:L){
  y = rnorm(1,Series[i-1],sigma)
  if(Uniform[i]<=exp(-abs(y)+abs(Series[i-1]))){
    Series[i] = y
    k=k+1
  }
  else{
    Series[i] = Series[i-1]
  }
}
AC_rate = k/L
print(AC_rate)
DATA_0 = as.data.frame(cbind(Series,1:L))
LP = function(x){1/2*exp(-abs(x))}
ggplot(data = DATA_0)+geom_histogram(mapping = aes(Series,..density..),bins = 70,fill='#FFCC00')+stat_function(fun=LP)+labs(title = "sigma=10")+theme(plot.title=element_text(hjust=0.5))+annotate(geom = 'text',x=5,y=0.4,label=as.character(AC_rate))
DATA_last = as.data.frame(cbind(Series[9000:10000],9000:10000))
ggplot(data = DATA_0)+geom_line(aes(V2,Series))+labs(title = "all steps")
ggplot(data=DATA_last)+geom_line(aes(V2,V1))+labs(title = "last 1000 steps")

## ----warning=FALSE,message=FALSE-----------------------------------------
library(ggplot2)
sigma = 1
L = 10000
Series = numeric(L)
Uniform = runif(L)
k=0
for (i in 2:L){
  y = rnorm(1,Series[i-1],sigma)
  if(Uniform[i]<=exp(-abs(y)+abs(Series[i-1]))){
    Series[i] = y
    k=k+1
  }
  else{
    Series[i] = Series[i-1]
  }
}
AC_rate = k/L
print(AC_rate)
DATA_0 = as.data.frame(cbind(Series,1:L))
LP = function(x){1/2*exp(-abs(x))}
ggplot(data = DATA_0)+geom_histogram(mapping = aes(Series,..density..),bins = 70,fill='#FFCC00')+stat_function(fun=LP)+labs(title = "sigma=1")+theme(plot.title=element_text(hjust=0.5))+annotate(geom = 'text',x=5,y=0.4,label=as.character(AC_rate))
DATA_last = as.data.frame(cbind(Series[9000:10000],9000:10000))
ggplot(data = DATA_0)+geom_line(aes(V2,Series))+labs(title = "all steps")
ggplot(data=DATA_last)+geom_line(aes(V2,V1))+labs(title = "last 1000 steps")

## ----warning=FALSE,message=FALSE-----------------------------------------
library(ggplot2)
sigma = 2
L = 10000
Series = numeric(L)
Uniform = runif(L)
k=0
for (i in 2:L){
  y = rnorm(1,Series[i-1],sigma)
  if(Uniform[i]<=exp(-abs(y)+abs(Series[i-1]))){
    Series[i] = y
    k=k+1
  }
  else{
    Series[i] = Series[i-1]
  }
}
AC_rate = k/L
print(AC_rate)
DATA_0 = as.data.frame(cbind(Series,1:L))
LP = function(x){1/2*exp(-abs(x))}
ggplot(data = DATA_0)+geom_histogram(mapping = aes(Series,..density..),bins = 70,fill='#FFCC00')+stat_function(fun=LP)+labs(title = "sigma=2")+theme(plot.title=element_text(hjust=0.5))+annotate(geom = 'text',x=5,y=0.4,label=as.character(AC_rate))
DATA_last = as.data.frame(cbind(Series[9000:10000],9000:10000))
ggplot(data = DATA_0)+geom_line(aes(V2,Series))+labs(title = "all steps")
ggplot(data=DATA_last)+geom_line(aes(V2,V1))+labs(title = "last 1000 steps")

## ----warning=FALSE,message=FALSE-----------------------------------------
library(ggplot2)
sigma = 5
L = 10000
Series = numeric(L)
Uniform = runif(L)
k=0
for (i in 2:L){
  y = rnorm(1,Series[i-1],sigma)
  if(Uniform[i]<=exp(-abs(y)+abs(Series[i-1]))){
    Series[i] = y
    k=k+1
  }
  else{
    Series[i] = Series[i-1]
  }
}
AC_rate = k/L
print(AC_rate)
DATA_0 = as.data.frame(cbind(Series,1:L))
LP = function(x){1/2*exp(-abs(x))}
ggplot(data = DATA_0)+geom_histogram(mapping = aes(Series,..density..),bins = 70,fill='#FFCC00')+stat_function(fun=LP)+labs(title = "sigma=5")+theme(plot.title=element_text(hjust=0.5))+annotate(geom = 'text',x=5,y=0.4,label=as.character(AC_rate))
DATA_last = as.data.frame(cbind(Series[9000:10000],9000:10000))
ggplot(data = DATA_0)+geom_line(aes(V2,Series))+labs(title = "all steps")
ggplot(data=DATA_last)+geom_line(aes(V2,V1))+labs(title = "last 1000 steps")

## ----warning=FALSE,message=FALSE-----------------------------------------
library(ggplot2)
sigma = 10
L = 10000
Series = numeric(L)
Uniform = runif(L)
k=0
for (i in 2:L){
  y = rnorm(1,Series[i-1],sigma)
  if(Uniform[i]<=exp(-abs(y)+abs(Series[i-1]))){
    Series[i] = y
    k=k+1
  }
  else{
    Series[i] = Series[i-1]
  }
}
AC_rate = k/L
print(AC_rate)
DATA_0 = as.data.frame(cbind(Series,1:L))
LP = function(x){1/2*exp(-abs(x))}
ggplot(data = DATA_0)+geom_histogram(mapping = aes(Series,..density..),bins = 70,fill='#FFCC00')+stat_function(fun=LP)+labs(title = "sigma=10")+theme(plot.title=element_text(hjust=0.5))+annotate(geom = 'text',x=5,y=0.4,label=as.character(AC_rate))
DATA_last = as.data.frame(cbind(Series[9000:10000],9000:10000))
ggplot(data = DATA_0)+geom_line(aes(V2,Series))+labs(title = "all steps")
ggplot(data=DATA_last)+geom_line(aes(V2,V1))+labs(title = "last 1000 steps")

## ----warning=FALSE,message=FALSE-----------------------------------------
data(mtcars)
formulas = list(mpg ~ disp,mpg ~ I(1 / disp),mpg ~ disp + wt,mpg ~ I(1 / disp) + wt)
# for loop
result_1 = list()
for (i in formulas){
  result_1 = c(result_1,list(lm(data=mtcars,i)))
}
#lapply
result_2 = lapply(formulas,lm,data=mtcars)
print(result_1)
print(result_2)

## ----warning=FALSE, message=FALSE----------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
#for loop
result1 = list()
for(i in bootstraps){
  result1 = c(result1,list(lm(data = i,mpg~disp)))
}
#lapply
result2 = lapply(bootstraps, lm,formula = mpg~disp)
print(result1)
print(result2)

#without anonymous function:
bootstraps = lapply(rep(list(mtcars),10),FUN = '[', i = sample(1:nrow(mtcars),rep=TRUE),j=1:ncol(mtcars) )
#for loop
result3 = list()
for(i in bootstraps){
  result3 = c(result3,list(lm(data = i,mpg~disp)))
}
#lapply
result4 = lapply(bootstraps, lm,formula = mpg~disp)
print(result3)
print(result4)

## ----warning=FALSE,message=FALSE-----------------------------------------
rsq <- function(mod) summary(mod)$r.squared
#question 1
R_1 = lapply(result_1,rsq)
R_2 = lapply(result_2,rsq)
#question 2
R1 = lapply(result1, rsq)
R2 = lapply(result2, rsq)
print(R_1)
print(R_2)
print(R1)
print(R2)

## ----warning=FALSE, message=FALSE----------------------------------------
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)
# Using anonymous function 
result_anonymous = sapply(trials, function(i) i$p.value)
# Don't use anonymous fnction
result_NO_anonymous = sapply(trials,'[[',i=3)
print(result_anonymous)
print(result_NO_anonymous)

## ----warning=FALSE, message=FALSE----------------------------------------
#mcsapply
mcsapply = function(X,FUN,mc.preschedule = TRUE, mc.set.seed = TRUE,mc.silent = FALSE, mc.cores = 1L,mc.cleanup = TRUE, mc.allow.recursive = TRUE,affinity.list = NULL,...){
  temp = parallel::mclapply(X=X,FUN=FUN,mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,mc.silent = mc.silent, mc.cores = mc.cores,mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive,affinity.list = affinity.list,... = ...)
  return(unlist(temp))
}
#mcvapply
mcvapply = function(X,FUN,mc.preschedule = TRUE, mc.set.seed = TRUE,mc.silent = FALSE, mc.cores = 1L,mc.cleanup = TRUE, mc.allow.recursive = TRUE,affinity.list = NULL,FUN.VALUE=NULL,...){
  temp = parallel::mclapply(X=X,FUN=FUN,mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,mc.silent = mc.silent, mc.cores = mc.cores,mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive,affinity.list = affinity.list,... = ...)
  is_FUN_VALUE = sapply(temp, function(i){
    a = (mode(i)==mode(FUN.VALUE))
    b = (length(i)==length(FUN.VALUE))
    return(a&b)
  })
  if(sum(is_FUN_VALUE)==length(is_FUN_VALUE)){
    return(unlist(temp))
  }
  else{
    stop("FUN.VALUE doesn't match the shape of return values")
  }
}

## ----warning=FALSE, message=FALSE----------------------------------------
library(microbenchmark)
library(Rcpp)
RF = function(L,sigma){
Series = numeric(L)
Uniform = runif(L)
k=0
for (i in 2:L){
  y = rnorm(1,Series[i-1],sigma)
  if(Uniform[i]<=exp(-abs(y)+abs(Series[i-1]))){
    Series[i] = y
    k=k+1
  }
  else{
    Series[i] = Series[i-1]
  }
}
return(Series)
}
#sourceCpp('CF.cpp')
#qqplot(RF(10000,1),CF(10000,1))
#COMP = summary(microbenchmark(RF(1000,1),CF(1000,1)))
#knitr::kable(COMP)

## ------------------------------------------------------------------------
# #Cpp function
# include<Rcpp.h>
# include<Rmath.h>
# include<cmath>
# using namespace Rcpp;
# //[[Rcpp::export]]
# NumericVector CF(int L, double sigma) {
#   NumericVector Series (L);
#   NumericVector Uniform = runif(L);
#   int k=0;
#   for (int i =1; i<L; i++){
#     NumericVector y = rnorm(1,Series[i-1],sigma);
#     if(Uniform[i]<=exp(-fabs(y[0])+fabs(Series[i-1]))){
#       Series[i] = y[0];
#       k=k+1;
#     }
#     else{
#       Series[i] = Series[i-1];
#     }
#   }
#   return Series;
# }

