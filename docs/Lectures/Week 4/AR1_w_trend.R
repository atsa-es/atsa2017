library(MARSS)
library(forecast)
alpha = 1
beta = .5
t = 1:100
trend = alpha + beta*t
#Trend with AR-1 errors
q = .02
b = .8
x = arima.sim(n=100, list(order=c(1,0,0), ar=b, sd=sqrt(q)))

#
mod.list=list(
  B=matrix("b"), U=matrix(0), Q=matrix("q"),
  Z=matrix(1), A=matrix("mean"), R=matrix(0),
  x0=matrix(x[1]), tinitx=1 )


fit=MARSS(as.vector(x), model=mod.list, method="BFGS")

fit.arima=arima(x, order=c(1,0,0))

mod.list=list(
  B=matrix(coef(fit.arima)[1]), U=matrix(0), Q=matrix(fit.arima$sigma2),
  Z=matrix(1), A=matrix(coef(fit.arima)[2]), R=matrix(0),
  x0=matrix(x[1]), tinitx=1 )

fit=MARSS(as.vector(x), model=mod.list)
