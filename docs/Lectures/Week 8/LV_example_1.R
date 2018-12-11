#------------------
set.seed(102)
#------------------

#Lotka-Volterra simulation
#(H_t+1-H_t) = b*dt*H(1-H/K)-a*dt*P*H 
#  = b*dt*H_t - b/K*dt*H_t*H_t- a*dt*P_t*H_t
# dP/dt = e(PaH)-sP
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
n=100;
x = matrix(0,2,simlen);
x[,1]=c(H_i, P_i);
rownames(x)=c("herbivore","predator")

for(i in 2:simlen){
x[1,i] = (b*dt*(1-x[1,i-1]/K)+rnorm(1,0,sd_H))*x[1,i-1] - a*dt*x[1,i-1]*x[2,i-1]+x[1,i-1];
x[2,i] = a*e*dt*x[1,i-1]*x[2,i-1]-s*dt*exp(rnorm(1,0,sd_P))*x[2,i-1]+x[2,i-1];
}

x=x[,(simlen-n+1):simlen];
x=log(x)

z.x = x-apply(x,1,mean)

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

LV.B = coef(LVfit, type="matrix")$B
rownames(LV.B)=rownames(x)
colnames(LV.B)=rownames(x)
abline(h=K, col="red", lty=2)

par(mfrow=c(2,1))
matplot(t(exp(x)), type="l", col=c("red","black"))
title(paste("red is herbivore; black is predator\nherb b[1,1]=",format(coef(LVfit)$B[1],digits=2)," conv eff=",e))
abline(h=K, col="red", lty=2)

#Lotka-Volterra simulation
#(H_t+1-H_t) = b*dt*H(1-H/K)-a*dt*P*H 
#  = b*dt*H_t - b/K*dt*H_t*H_t- a*dt*P_t*H_t
# dP/dt = e(PaH)-sP
a=0.7; #attack rate
e=0.8; #conv eff
K=5;
dt=.2;
b=4; #birth rate of herbivore
s=0.5; #death rate of predator
H_i=5; P_i=7;
sd_H = sqrt(.02);
sd_P = sqrt(.02);

simlen = 1000;
n=100;
x = matrix(0,2,simlen);
x[,1]=c(H_i, P_i);
rownames(x)=c("herbivore","predator")

for(i in 2:simlen){
  x[1,i] = (b*dt*(1-x[1,i-1]/K)+rnorm(1,0,sd_H))*x[1,i-1] - a*dt*x[1,i-1]*x[2,i-1]+x[1,i-1];
  x[2,i] = a*e*dt*x[1,i-1]*x[2,i-1]-s*dt*exp(rnorm(1,0,sd_P))*x[2,i-1]+x[2,i-1];
}

x=x[,(simlen-n+1):simlen];
x=log(x)

z.x = x-apply(x,1,mean)

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

LV.B = coef(LVfit, type="matrix")$B
rownames(LV.B)=rownames(x)
colnames(LV.B)=rownames(x)

matplot(t(exp(x)), type="l", col=c("red","black"))
title(paste("red is herbivore; black is predator\nherb b[1,1]=",format(coef(LVfit)$B[1],digits=2)," conv eff=",e))
abline(h=K, col="red", lty=2)
