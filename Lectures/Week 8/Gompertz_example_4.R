#################################################################
##  Discrete Gompertz
##  Example 4.  Comparison of ML and REML estimates
#################################################################

##A mean-reverting model
K=log(2000)
u=1
b=1-u/K
x0=K
q=.02 #increased

#generate and underlying true process
x=b*x0+u+rnorm(1,0,sqrt(q))
n=40
for(i in 2:n) x[i]=b*x[i-1]+u+rnorm(1,0,sqrt(q))

#add some observation error
bias = -20
r=.1
y = x + rnorm(n,0,sqrt(r)) + bias

#Dennis et al 2006
Observed.t=exp(y)
a0=u # Initial value of a 
c0=b # Initial value of c 
ssq0=q # Initial value of ssq 
tsq0=r # Initial value of tsq
source("GompertzSS_REML.R")

#MARSS unconstrained
mod.list=list(
  U=matrix("u"),
  x0=matrix("x0"),
  B=matrix("b"),
  Q=matrix("q"),
  Z=matrix(1),
  A=matrix(0),
  R=matrix("r"),
  tinitx=1)

#Notice that R goes to 0 in this example and b is small
fit1=MARSS(y,model=mod.list,method="BFGS")

#MARSS constrained
#De-mean and set U=0
mod.list=list(
  U=matrix(0),
  x0=matrix("x0"),
  B=matrix("b"),
  Q=matrix("q"),
  Z=matrix(1),
  A=matrix(0),
  R=matrix("r"),
  tinitx=1)

fit2=MARSS(y-mean(y),model=mod.list,method="BFGS")
#fit2b=MARSS(y-mean(y),model=mod.list)

# Example using redstart data from Dennis et al 2006
# This is an estimate where the time series is not assumed to be stationary
redstart=c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2, 4,11,11,9,6) 
y=log(redstart)
plot(y, main="redstart counts")
#BFGS gets hung up on edge!
fit.red=MARSS(y-mean(y),model=mod.list,method="BFGS")
fit.red=MARSS(y-mean(y),model=mod.list)

#Dennis et al 2006
Observed.t=exp(y)
a0=u # Initial value of a 
c0=b # Initial value of c 
ssq0=q # Initial value of ssq 
tsq0=r # Initial value of tsq
source("GompertzSS_REML.R")

#But!
#Dennis et al 2006
Observed.t=exp(y+10)
a0=u # Initial value of a 
c0=b # Initial value of c 
ssq0=q # Initial value of ssq 
tsq0=r # Initial value of tsq
source("GompertzSS_REML.R")

