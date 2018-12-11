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
H_i=K; P_i=7;
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

#Add error
y=x+t(mvrnorm(n,rep(0,2),diag(.1,2)))

z.y = y - apply(y,1,mean)
mod.list=list(
  B="unconstrained",
  U="zero",
  Q="diagonal and unequal",
  Z="identity",
  A="zero",
  R="diagonal and equal",
  x0="unequal",
  tinitx=1)

LVfit.werr = MARSS(z.y,model=mod.list, method="BFGS")

LV.B = coef(LVfit.werr, type="matrix")$B
rownames(LV.B)=rownames(x)
colnames(LV.B)=rownames(x)
print(LV.B)
