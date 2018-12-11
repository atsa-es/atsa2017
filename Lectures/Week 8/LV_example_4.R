#------------------

#set.seed(102)
#------------------

#Lotka-Volterra simulation
# Covariate affects a (attack rate)
################################################
# Sinusoidal covariate drives K
################################################

#(H_t+1-H_t) = b*dt*H(1-H/K)-a*dt*P*H 
#  = b*dt*H_t - b/K*dt*H_t*H_t- a*dt*P_t*H_t
# dP/dt = e(PaH)-sP
n=100;
a=0.7; #attack rate
e0=0.7; #conv eff
K=5;
dt=.2;
b=4; #birth rate of herbivore
s=0.5; #death rate of predator
H_i=K; P_i=7;
sd_H = sqrt(.02);
sd_P = sqrt(.02);

simlen = 1000;

covariate=.25*sin(.1*pi*1:simlen)
e=e0*exp(covariate)

x = matrix(NA,2,simlen);
x[,1]=c(H_i, P_i);
rownames(x)=c("herbivore","predator")

#While so that reruns if any NAs in sim data
while(any(is.na(x)) | any(is.nan(x))){
  for(i in 2:simlen){
    x[1,i] = (b*dt*(1-x[1,i-1]/K)+rnorm(1,0,sd_H))*x[1,i-1] - a*dt*x[1,i-1]*x[2,i-1]+x[1,i-1];
    x[2,i] = a*e[i-1]*dt*x[1,i-1]*x[2,i-1]-s*dt*exp(rnorm(1,0,sd_P))*x[2,i-1]+x[2,i-1];
  }
}

x=x[,(simlen-n+1):simlen];
x=log(x)
c=covariate[(simlen-n+1):simlen]

#demean
z.x = (x-apply(x,1,mean))/sqrt(apply(x,1,var))
z.c = (c - mean(c))/sqrt(var(c))

## Print the output
par(mfrow=c(2,1))
plot(c,type="l",main="covariate")
ylims=c(0,max(max(K),max(exp(x))))
matplot(t(exp(x)), type="l", col=c("red","black"),ylim=ylims,lty=1)
abline(h=K, col="red", lty=2)


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

title("red is herbivore; black is predator")

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
cat("\nSinusoidal covariate affects conversion efficiency of predator; fit with covariate\n\n");print(LV.B2);cat("\nfit WITHOUT covariate\n\n");print(LV.B)

readline("Continue?")


################################################
# Linear covariate drives e
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
e0=0.7; #conv eff
K=5;
dt=.2;
b=4; #birth rate of herbivore
s=0.5; #death rate of predator
H_i=5; P_i=7;
sd_H = sqrt(.02);
sd_P = sqrt(.02);

simlen = 1000;

covariate=c(rep(1,simlen-n),1-.01*1:n)
e=e0*covariate

x = matrix(NA,2,simlen);
x[,1]=c(H_i, P_i);
rownames(x)=c("herbivore","predator")

#While so that reruns if any NAs in sim data
while(any(is.na(x)) | any(is.nan(x))){
  for(i in 2:simlen){
  x[1,i] = (b*dt*(1-x[1,i-1]/K)+rnorm(1,0,sd_H))*x[1,i-1] - a*dt*x[1,i-1]*x[2,i-1]+x[1,i-1];
  x[2,i] = a*e[i-1]*dt*x[1,i-1]*x[2,i-1]-s*dt*exp(rnorm(1,0,sd_P))*x[2,i-1]+x[2,i-1];
}
}

x=x[,(simlen-n+1):simlen];
x=log(x)
c=covariate[(simlen-n+1):simlen]

#demean
z.x = (x-apply(x,1,mean))/sqrt(apply(x,1,var))
z.c = (c - mean(c))/sqrt(var(c))

## Print the output
par(mfrow=c(2,1))
plot(c,type="l",main="covariate")

ylims=c(0,max(max(K),max(exp(x))))
matplot(t(exp(x)), type="l", col=c("red","black"),ylim=ylims,lty=1)
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

title("red is herbivore; black is predator")

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
cat("\nLinear covariate affects e of predator; fit with covariate\n\n");print(LV.B2);cat("\nfit WITHOUT covariate\n\n");print(LV.B)

readline("Continue?")


################################################
# Random Walk with drift drives e
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
e0=0.2; #conv eff
K=5;
dt=.2;
b=4; #birth rate of herbivore
s=0.5; #death rate of predator
H_i=5; P_i=7;
sd_H = sqrt(.02);
sd_P = sqrt(.02);

simlen = 1000;

covariate=c(rep(0,simlen-n),cumsum(rnorm(n,-.02,.1)))
covariate=1+.5*covariate/sqrt(var(covariate))
e=e0*exp(covariate)

x = matrix(NA,2,simlen);
x[,1]=c(H_i, P_i);
rownames(x)=c("herbivore","predator")

#While so that reruns if any NAs in sim data
while(any(is.na(x)) | any(is.nan(x))){
  for(i in 2:simlen){
    x[1,i] = (b*dt*(1-x[1,i-1]/K)+rnorm(1,0,sd_H))*x[1,i-1] - a*dt*x[1,i-1]*x[2,i-1]+x[1,i-1];
    x[2,i] = a*e[i-1]*dt*x[1,i-1]*x[2,i-1]-s*dt*exp(rnorm(1,0,sd_P))*x[2,i-1]+x[2,i-1];
  }
}

x=x[,(simlen-n+1):simlen];
x=log(x)
c=covariate[(simlen-n+1):simlen]

#demean
z.x = (x-apply(x,1,mean))/sqrt(apply(x,1,var))
z.c = (c - mean(c))/sqrt(var(c))

## Print the output
par(mfrow=c(2,1))
plot(c,type="l",main="covariate")

ylims=c(0,max(max(K),max(exp(x))))
matplot(t(exp(x)), type="l", col=c("red","black"),ylim=ylims,lty=1)
abline(h=K, col="red", lty=2)

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

title("red is herbivore; black is predator")

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
cat("\nRW covariate affects e of predator; fit with covariate\n\n");print(LV.B2);cat("\nfit WITHOUT covariate\n\n");print(LV.B)

readline("Continue?")


################################################
# Mean-reverting random Noise drives e
################################################

#(H_t+1-H_t) = b*dt*H(1-H/K)-a*dt*P*H 
#  = b*dt*H_t - b/K*dt*H_t*H_t- a*dt*P_t*H_t
# dP/dt = e(PaH)-sP
n=100;
a=0.7; #attack rate
e0=0.2; #conv eff
K=5;
dt=.2;
b=4; #birth rate of herbivore
s=0.5; #death rate of predator
H_i=5; P_i=7;
sd_H = sqrt(.02);
sd_P = sqrt(.02);

simlen = 1000;

covariate=arima.sim(simlen,model=list(ar=.8,sd=.05))
covariate=1+.5*covariate/sqrt(var(covariate))
e=e0*exp(covariate)

x = matrix(NA,2,simlen);
x[,1]=c(H_i, P_i);
rownames(x)=c("herbivore","predator")

#While so that reruns if any NAs in sim data
while(any(is.na(x)) | any(is.nan(x))){
  for(i in 2:simlen){
    x[1,i] = (b*dt*(1-x[1,i-1]/K)+rnorm(1,0,sd_H))*x[1,i-1] - a*dt*x[1,i-1]*x[2,i-1]+x[1,i-1];
    x[2,i] = a*e[i-1]*dt*x[1,i-1]*x[2,i-1]-s*dt*exp(rnorm(1,0,sd_P))*x[2,i-1]+x[2,i-1];
  }
}

x=x[,(simlen-n+1):simlen];
x=log(x)
c=covariate[(simlen-n+1):simlen]

#demean and standardize
z.x = (x-apply(x,1,mean))/sqrt(apply(x,1,var))
z.c = (c - mean(c))/sqrt(var(c))

## Print the output
par(mfrow=c(2,1))
plot(c,type="l",main="covariate")

ylims=c(0,max(max(K),max(exp(x))))
matplot(t(exp(x)), type="l", col=c("red","black"),ylim=ylims,lty=1)
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

title("red is herbivore; black is predator")

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
cat("\nMean-reverting random Walk noise covariate affects e of predator; fit with covariate\n\n");print(LV.B2);cat("\nfit WITHOUT covariate\n\n");print(LV.B)


