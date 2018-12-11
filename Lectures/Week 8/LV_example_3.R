#Lotka-Volterra simulation
# Covariate affects K (herbivore K)
# You MUST demean the covariates! Just like you demean data, otherwise setting u=zero doesn't work (mathematically).

#------------------
#set.seed(102)
#------------------

################################################
# Sinusoidal covariate drives K
################################################

#(H_t+1-H_t) = b*dt*H(1-H/K)-a*dt*P*H 
#  = b*dt*H_t - b/K*dt*H_t*H_t- a*dt*P_t*H_t
# dP/dt = e(PaH)-sP
n=100;
a=0.7; #attack rate
e=0.2; #conv eff
K=5;
dt=.2;
b=4; #birth rate of herbivore
s=0.5; #death rate of predator
H_i=K; P_i=7;
sd_H = sqrt(.02);
sd_P = sqrt(.02);

simlen = 1000;

covariate=.25*sin(.1*pi*1:simlen)
K=K*exp(covariate)

x = matrix(NA,2,simlen);
x[,1]=c(H_i, P_i);
rownames(x)=c("herbivore","predator")

while(any(is.nan(x)) | any(is.na(x))){
for(i in 2:simlen){
x[1,i] = (b*dt*(1-x[1,i-1]/K[i-1])+rnorm(1,0,sd_H))*x[1,i-1] - a*dt*x[1,i-1]*x[2,i-1]+x[1,i-1];
x[2,i] = a*e*dt*x[1,i-1]*x[2,i-1]-s*dt*exp(rnorm(1,0,sd_P))*x[2,i-1]+x[2,i-1];
}
}

x=x[,(simlen-n+1):simlen];
x=log(x)
c=covariate[(simlen-n+1):simlen]

#demean data and covariate
z.x = x-apply(x,1,mean)
z.c = c-mean(c)

## Print the output
par(mfrow=c(2,1))
plot(c,type="l",main="covariate")

ylims=c(0,max(max(K),max(exp(x))))
matplot(t(exp(x)), type="l", col=c("red","black"),ylim=ylims,lty=1)
title("red is herbivore; black is predator")
lines(K[(simlen-n+1):simlen], col="red", lty=2)

## Ignore the covariate
mod.list=list(
  B="unconstrained",
  U="zero",
  Q="diagonal and unequal",
  Z="identity",
  A="zero",
  R="zero",
  x0=z.x[,1,drop=FALSE],
  tinitx=1)

LVfit = MARSS(z.x,model=mod.list)

## Include the covariate
mod.list2=list(
  B="unconstrained",
  U="zero",
  Q="diagonal and unequal",
  Z="identity",
  A="zero",
  R="zero",
  C="unconstrained",
  c=matrix(z.c,1,n),
  x0=z.x[,1,drop=FALSE],
  tinitx=1)

LVfit2 = MARSS(z.x,model=mod.list2)

#compare the 2 approaches
LV.B2 = cbind(coef(LVfit2, type="matrix")$B,coef(LVfit2,type="matrix")$C)
rownames(LV.B2)=rownames(x)
colnames(LV.B2)=c(rownames(x),"cov")

LV.B = coef(LVfit, type="matrix")$B
rownames(LV.B)=rownames(x)
colnames(LV.B)=rownames(x)
cat("\nSinusoidal covariate affects K of herbivore; fit with covariate\n\n");print(LV.B2);cat("\nfit WITHOUT covariate\n\n");print(LV.B)

readline("Continue?")

################################################
# Linear covariate drives K
# It took me awhile to find working param values
# a) Don't drive the predator extinct by having K too low
# b) Have K start big enough that DD is low then end low
#    enought that DD is big
################################################

#(H_t+1-H_t) = b*dt*H(1-H/K)-a*dt*P*H 
#  = b*dt*H_t - b/K*dt*H_t*H_t- a*dt*P_t*H_t
# dP/dt = e(PaH)-sP
n=100;
a=0.7; #attack rate
e=0.2; #conv eff
K=5;
dt=.2;
b=4; #birth rate of herbivore
s=0.5; #death rate of predator
H_i=5; P_i=7;
sd_H = sqrt(.02);
sd_P = sqrt(.02);

simlen = 1000;

covariate=c(rep(1,simlen-n),1-.0025*1:n)
K=K*covariate

x = matrix(NA,2,simlen);
x[,1]=c(H_i, P_i);
rownames(x)=c("herbivore","predator")

while(any(is.nan(x)) | any(is.na(x))){
  for(i in 2:simlen){
    x[1,i] = (b*dt*(1-x[1,i-1]/K[i-1])+rnorm(1,0,sd_H))*x[1,i-1] - a*dt*x[1,i-1]*x[2,i-1]+x[1,i-1];
    x[2,i] = a*e*dt*x[1,i-1]*x[2,i-1]-s*dt*exp(rnorm(1,0,sd_P))*x[2,i-1]+x[2,i-1];
  }
}

x=x[,(simlen-n+1):simlen];
x=log(x)
c=covariate[(simlen-n+1):simlen]

#demean data and covariate
z.x = x-apply(x,1,mean)
z.c = c-mean(c)

## Print the output
par(mfrow=c(2,1))
plot(c,type="l",main="covariate")

ylims=c(0,max(max(K),max(exp(x))))
matplot(t(exp(x)), type="l", col=c("red","black"),ylim=ylims,lty=1)
lines(K[(simlen-n+1):simlen], col="red", lty=2)
title("red is herbivore; black is predator")

## Ignore the covariate
mod.list=list(
  B="unconstrained",
  U="zero",
  Q="diagonal and unequal",
  Z="identity",
  A="zero",
  R="zero",
  x0=z.x[,1,drop=FALSE],
  tinitx=1)

LVfit = MARSS(z.x,model=mod.list)

## Include the covariate
mod.list2=list(
  B="unconstrained",
  U="zero",
  Q="diagonal and unequal",
  Z="identity",
  A="zero",
  R="zero",
  C="unconstrained",
  c=matrix(z.c,1,n),
  x0=z.x[,1,drop=FALSE],
  tinitx=1)

LVfit2 = MARSS(z.x,model=mod.list2)

#compare the 2 approaches
LV.B2 = cbind(coef(LVfit2, type="matrix")$B,coef(LVfit2,type="matrix")$C)
rownames(LV.B2)=rownames(x)
colnames(LV.B2)=c(rownames(x),"cov")

LV.B = coef(LVfit, type="matrix")$B
rownames(LV.B)=rownames(x)
colnames(LV.B)=rownames(x)
cat("\nLinear covariate affects K of herbivore; fit with covariate\n\n");print(LV.B2);cat("\nfit WITHOUT covariate\n\n");print(LV.B)

readline("Continue?")

################################################
# Random Walk with drift drives K
# It took me awhile to find working param values
# a) Don't drive the predator extinct by having K too low
# b) Have K start big enough that DD is low then end low
#    enought that DD is big
################################################

#(H_t+1-H_t) = b*dt*H(1-H/K)-a*dt*P*H 
#  = b*dt*H_t - b/K*dt*H_t*H_t- a*dt*P_t*H_t
# dP/dt = e(PaH)-sP
n=100;
a=0.7; #attack rate
e=0.2; #conv eff
K0=5;
dt=.2;
b=4; #birth rate of herbivore
s=0.5; #death rate of predator
H_i=5; P_i=7;
sd_H = sqrt(.02);
sd_P = sqrt(.02);

simlen = 1000;

covariate=c(rep(1,simlen-n),1+cumsum(rnorm(n,-.02,.1)))
K=K0*exp(covariate)

x = matrix(NA,2,simlen);
x[,1]=c(H_i, P_i);
rownames(x)=c("herbivore","predator")

while(any(is.nan(x)) | any(is.na(x))){
  for(i in 2:simlen){
    x[1,i] = (b*dt*(1-x[1,i-1]/K[i-1])+rnorm(1,0,sd_H))*x[1,i-1] - a*dt*x[1,i-1]*x[2,i-1]+x[1,i-1];
    x[2,i] = a*e*dt*x[1,i-1]*x[2,i-1]-s*dt*exp(rnorm(1,0,sd_P))*x[2,i-1]+x[2,i-1];
  }
}

x=x[,(simlen-n+1):simlen];
x=log(x)
c=covariate[(simlen-n+1):simlen]

#demean
z.x = x-apply(x,1,mean)
z.c = c-mean(c)

## Print the output
par(mfrow=c(2,1))
plot(c,type="l",main="covariate")

ylims=c(0,max(max(K[(simlen-n+1):simlen]),max(exp(x))))
matplot(t(exp(x)), type="l", col=c("red","black"),ylim=ylims,lty=1)
lines(K[(simlen-n+1):simlen], col="red", lty=2)
title("red is herbivore; black is predator")

## Ignore the covariate
mod.list=list(
  B="unconstrained",
  U="zero",
  Q="diagonal and unequal",
  Z="identity",
  A="zero",
  R="zero",
  x0=z.x[,1,drop=FALSE],
  tinitx=1)

LVfit = MARSS(z.x,model=mod.list)

## Include the covariate
mod.list2=list(
  B="unconstrained",
  U="zero",
  Q="diagonal and unequal",
  Z="identity",
  A="zero",
  R="zero",
  C="unconstrained",
  c=matrix(z.c,1,n),
  x0=z.x[,1,drop=FALSE],
  tinitx=1)

LVfit2 = MARSS(z.x,model=mod.list2)

#compare the 2 approaches
LV.B2 = cbind(coef(LVfit2, type="matrix")$B,coef(LVfit2)$C)
rownames(LV.B2)=rownames(x)
colnames(LV.B2)=c(rownames(x),"cov")

LV.B = coef(LVfit, type="matrix")$B
rownames(LV.B)=rownames(x)
colnames(LV.B)=rownames(x)
cat("\nRandom Walk covariate affects K of herbivore; fit with covariate\n\n");print(LV.B2);cat("\nfit WITHOUT covariate\n\n");print(LV.B)

readline("Continue?")

################################################
# Mean-reverting Random Noise drives K
################################################

#(H_t+1-H_t) = b*dt*H(1-H/K)-a*dt*P*H 
#  = b*dt*H_t - b/K*dt*H_t*H_t- a*dt*P_t*H_t
# dP/dt = e(PaH)-sP
n=100;
a=0.7; #attack rate
e=0.2; #conv eff
K0=5;
dt=.2;
b=4; #birth rate of herbivore
s=0.5; #death rate of predator
H_i=5; P_i=7;
sd_H = sqrt(.02);
sd_P = sqrt(.02);

simlen = 1000;

covariate=arima.sim(simlen,model=list(ar=.8,sd=.05))
covariate=1+.5*covariate/sqrt(var(covariate))
K=K0*exp(covariate)

x = matrix(NA,2,simlen);
x[,1]=c(H_i, P_i);
rownames(x)=c("herbivore","predator")

while(any(is.nan(x)) | any(is.na(x))){
  for(i in 2:simlen){
    x[1,i] = (b*dt*(1-x[1,i-1]/K[i-1])+rnorm(1,0,sd_H))*x[1,i-1] - a*dt*x[1,i-1]*x[2,i-1]+x[1,i-1];
    x[2,i] = a*e*dt*x[1,i-1]*x[2,i-1]-s*dt*exp(rnorm(1,0,sd_P))*x[2,i-1]+x[2,i-1];
  }
}

x=x[,(simlen-n+1):simlen];
x=log(x)
c=covariate[(simlen-n+1):simlen]

#demean
z.x = x-apply(x,1,mean)
z.c = c - mean(c)

## Print the output
par(mfrow=c(2,1))
plot(c,type="l",main="covariate")

ylims=c(0,max(max(K[(simlen-n+1):simlen]),max(exp(x))))
matplot(t(exp(x)), type="l", col=c("red","black"),ylim=ylims,lty=1)
lines(K[(simlen-n+1):simlen], col="red", lty=2)
title("red is herbivore; black is predator")

## Ignore the covariate
mod.list=list(
  B="unconstrained",
  U="zero",
  Q="diagonal and unequal",
  Z="identity",
  A="zero",
  R="zero",
  x0=z.x[,1,drop=FALSE],
  tinitx=1)

LVfit = MARSS(z.x,model=mod.list)

## Include the covariate
mod.list2=list(
  B="unconstrained",
  U="zero",
  Q="diagonal and unequal",
  Z="identity",
  A="zero",
  R="zero",
  C="unconstrained",
  c=matrix(z.c,1,n),
  x0=z.x[,1,drop=FALSE],
  tinitx=1)

LVfit2 = MARSS(z.x,model=mod.list2)

#compare the 2 approaches
LV.B2 = cbind(coef(LVfit2, type="matrix")$B,coef(LVfit2)$C)
rownames(LV.B2)=rownames(x)
colnames(LV.B2)=c(rownames(x),"cov")

LV.B = coef(LVfit, type="matrix")$B
rownames(LV.B)=rownames(x)
colnames(LV.B)=rownames(x)
cat("\nMean-reverting random walk covariate affects K of herbivore; fit with covariate\n\n");print(LV.B2);cat("\nfit WITHOUT covariate\n\n");print(LV.B)
