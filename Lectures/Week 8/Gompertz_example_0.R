#################################################################
##  Discrete Gompertz
##  Example 0.  Observation error leads to spurious dens-dep
#################################################################

require(MARSS)
library(forecast)

#-----
set.seed(123)
#-----

##A mean-reverting model
K=log(2000)
u=1
b=1-u/K
x0=1
q=.02

#generate an underlying true process
#x(1)
x=b*x0+u+rnorm(1,0,sqrt(q))
#add x(2:n)
n=40
for(i in 2:n) x[i]=b*x[i-1]+u+rnorm(1,0,sqrt(q))

#No observation error
y=x


#fit with Arima
fit.Arima=Arima(y, order=c(1,0,0))

# true versus estimated
c(true=b, est=coef(fit.Arima)["ar1"])

#that's not so good why did that happen?

#Plot y
plot(y, type="l")

#fit with MARSS
mod.list=list(
  U=matrix("u"),
  x0=matrix(y[1]),
  B=matrix("b"),
  Q=matrix("q"),
  Z=matrix(1),
  A=matrix(0),
  R=matrix(0),
  tinitx=1)

#Note u and b estimation is hard; you might need to use method="BFGS"
fit.marss=MARSS(y,model=mod.list, method="BFGS")

# true versus estimated
c(true=b, 
  est.Arima=coef(fit.Arima)["ar1"], 
  est.marss=coef(fit.marss)$B)

#Let's add some observation error
r=.1
y = x + rnorm(n,0,sqrt(r))


mod.list=list(
  U=matrix("u"),
  x0=matrix(y[1]),
  B=matrix("b"),
  Q=matrix("q"),
  Z=matrix(1),
  A=matrix(0),
  R=matrix(0),
  tinitx=1)

#fit
fit.Arima2=Arima(y, order=c(1,0,0))
fit.marss2=MARSS(y,model=mod.list, method="BFGS")

# true versus estimated
c(true=b, 
  est.Arima=coef(fit.Arima2)["ar1"], 
  est.marss=coef(fit.marss2)$B)

#Look at the acfs

#-----
par.old=par(mfrow=c(3,1), mai=c(0.4,0.7,0.1,0.1))
#-----

#take a look at the acf of the diff( white noise )
acf(x) #acf of AR-1 b!=1
acf(diff(x)) #acf of diff of that
acf(diff(rnorm(100))) #characteristic of diff of iid error

#take a look at the acf of the diff( white noise )
acf(y) #acf of AR-1 b!=1
acf(diff(y)) #acf of diff of that
acf(diff(rnorm(100))) #characteristic of diff of iid error

plot(1:n,exp(x),xlim=c(0,n),type="l",lwd=2)
points(exp(y))
acf(diff(x),ylab="ACF diff x")
acf(diff(y),ylab="ACF diff y")

par(par.old)
