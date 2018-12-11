### R code from vignette source 'fitting-univariate-ss-key-content.Rnw'

###################################################
### code chunk number 1: RUNFIRST
###################################################
require(MARSS)
options(prompt=" ", continue=" ", width=60)


###################################################
### code chunk number 2: key-set-up-data
###################################################
library(MARSS)
dat=log(grouse[,2])


###################################################
### code chunk number 3: hw1-fig
###################################################
plot(grouse[,1], dat, type="l", ylab="log count", xlab="")


###################################################
### code chunk number 4: hw1-fig-plot
###################################################
par(mar=c(2, 4, 2, 2))
plot(grouse[,1], dat, type="l", ylab="log count", xlab="")


###################################################
### code chunk number 5: mod-lists
###################################################
mod.list1=list(
  B=matrix(1), U=matrix(0), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  x0=matrix("a"), tinitx=0)
fit1.marss = MARSS(dat, model=mod.list1)

mod.list2=list(
  B=matrix(1), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  x0=matrix("a"), tinitx=0)
fit2.marss = MARSS(dat, model=mod.list2)


###################################################
### code chunk number 6: mod-res
###################################################
coef(fit1.marss, type="vector")
coef(fit2.marss, type="vector")


###################################################
### code chunk number 7: mod-aics
###################################################
c(fit1.marss$AICc, fit2.marss$AICc)


###################################################
### code chunk number 8: mod-auto-arima
###################################################
library(forecast)
auto.arima(dat)


###################################################
### code chunk number 9: mod-auto-arima-trace
###################################################
auto.arima(dat, trace=TRUE)


###################################################
### code chunk number 10: comp-arima-aicc
###################################################
fit1.arima=Arima(dat, order=c(0,1,0))
fit2.arima=Arima(dat, order=c(0,1,0), include.drift=TRUE)

fit2.arima$aicc-fit1.arima$aicc
fit2.marss$AICc-fit1.marss$AICc


###################################################
### code chunk number 11: hw2.data
###################################################
dat=cumsum(rnorm(100,0.1,1))


###################################################
### code chunk number 12: hw2-arima
###################################################
fit.arima=Arima(dat, order=c(0,1,0), include.drift=TRUE)


###################################################
### code chunk number 13: hw2-marss
###################################################
mod.list=list(
  B=matrix(1), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  x0=matrix("mu"), tinitx=0)
fit.marss = MARSS(dat, model=mod.list)


###################################################
### code chunk number 14: hw2-marss-alt
###################################################
mod.list.alt=list(
  B=matrix(1), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  x0=matrix(dat[1]), tinitx=1)
fit.alt.marss = MARSS(dat, model=mod.list.alt, method="BFGS")


###################################################
### code chunk number 15: hw2-compare-coef
###################################################
coef(fit.marss, type="vector")
coef(fit.alt.marss, type="vector")
c(coef(fit.arima), s2=fit.arima$sigma2)


###################################################
### code chunk number 16: hw3.data
###################################################
diff.dat=diff(dat)


###################################################
### code chunk number 17: hw2-fit-diff-dat
###################################################
fit.diff.arima=Arima(diff.dat, order=c(0,0,0), include.mean=TRUE)


###################################################
### code chunk number 18: hw2-fit-diff-marss
###################################################
mod.list.diff.1=list(
  B=matrix(0), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  x0=matrix(0), tinitx=0)
fit.alt.diff.1 = MARSS(diff.dat, model=mod.list.diff.1)



###################################################
### code chunk number 19: hw2-fit-diff-marss-2
###################################################
mod.list.diff.2=list(
  B=matrix(0), U=matrix(0), Q=matrix(0),
  Z=matrix(0), A=matrix("u"), R=matrix("r"),
  x0=matrix(0), tinitx=0)
fit.alt.diff.2 = MARSS(diff.dat, model=mod.list.diff.2)



###################################################
### code chunk number 20: hw2-fit-diff-coef
###################################################
coef(fit.alt.diff.1, type="vector")
coef(fit.alt.diff.2, type="vector")
c(coef(fit.diff.arima), s2=fit.diff.arima$sigma2)


###################################################
### code chunk number 21: hw3-sim
###################################################
#set up my parameter values
b=.8; u=2; x0=10; q=0.1
nsim=1000
#set up my holder for x
x=rep(NA, nsim)
x[1]=b*x0+u+rnorm(1,0,sqrt(q))
for(t in 2:nsim) x[t]=b*x[t-1]+u+rnorm(1,0,sqrt(q))


###################################################
### code chunk number 22: hw3-fig
###################################################
plot(x, type="l",xlab="", ylab="x")


###################################################
### code chunk number 23: hw3-fig-plot
###################################################
par(mar=c(2, 4, 2, 2))
plot(x, type="l",xlab="", ylab="x", ylim=c(-6+u/(1-b),6+u/(1-b)))


###################################################
### code chunk number 24: hw3-sim
###################################################
#set up my parameter values
u2=u+1
x2=rep(NA, nsim)
x2[1]=b*x0+u2+rnorm(1,0,sqrt(q))
for(t in 2:nsim) x2[t]=b*x2[t-1]+u2+rnorm(1,0,sqrt(q))
#second u
u3=u-1
x3=rep(NA, nsim)
x3[1]=b*x0+u3+rnorm(1,0,sqrt(q))
for(t in 2:nsim) x3[t]=b*x3[t-1]+u3+rnorm(1,0,sqrt(q))



###################################################
### code chunk number 25: hw3-fig-plot2
###################################################
par(mar=c(2, 4, 2, 2))
plot(x, type="l",xlab="", ylab="x", ylim=c(-6+u/(1-b),6+u/(1-b)))
lines(x2, col="blue")
lines(x3, col="red")
legend("bottomright", c("u+1","u-1"), col=c("blue", "red"), lty=1, bg="white")


###################################################
### code chunk number 26: hw3-sim2
###################################################
#set up my parameter values
b1=0.9
x0=u/(1-b1)
x1=rep(NA, nsim)
x1[1]=b1*x0+u+rnorm(1,0,sqrt(q))
for(t in 2:nsim) x1[t]=b1*x1[t-1]+u+rnorm(1,0,sqrt(q))
# second b
b2=0.1
x0=u/(1-b2)
x2=rep(NA, nsim)
x2[1]=b2*x0+u+rnorm(1,0,sqrt(q))
for(t in 2:nsim) x2[t]=b2*x2[t-1]+u+rnorm(1,0,sqrt(q))


###################################################
### code chunk number 27: hw3-fig-plot3
###################################################
par(mfrow=c(1,2), mar=c(2,2,2,1))
plot(x1, type="l",xlab="", ylab="x", main="b=0.9")
plot(x2, type="l",xlab="", ylab="x", main="b=0.1")


###################################################
### code chunk number 28: hw3-sim3
###################################################
#set up my parameter values
b=0.9
u=1
x0=u/(1-b)
err=rnorm(nsim,0,sqrt(q))
x1=rep(NA, nsim)
x1[1]=b*x0+u+err[1]
for(t in 2:nsim) x1[t]=b*x1[t-1]+u+err[t]
# second u
u=2
x0=u/(1-b)
x2=rep(NA, nsim)
x2[1]=b*x0+u+err[1]
for(t in 2:nsim) x2[t]=b*x2[t-1]+u+err[t]


###################################################
### code chunk number 29: hw3-fig-plot4
###################################################
par(mfrow=c(1,2), mar=c(2,2,2,1))
plot(x1, type="l",xlab="", ylab="x", main="u=1")
plot(x2, type="l",xlab="", ylab="x", main="u=2")


###################################################
### code chunk number 30: hw4
###################################################
library(MARSS)
dat=log(graywhales[,2])


###################################################
### code chunk number 31: hw4-fit
###################################################
mod.list=list(
  B=matrix(1), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix("r"),
  x0=matrix("mu"), tinitx=0)
fit.marss = MARSS(dat, model=mod.list)


###################################################
### code chunk number 32: hw4-fig-plot1
###################################################
par(mar=c(2,2,2,2))
plot(graywhales[,1], fit.marss$states[1,], type="l",xlab="", ylab="log count")
points(graywhales[,1], dat)


###################################################
### code chunk number 33: hw4-sim
###################################################
#1997 is the 39th (last) data point
x0=fit.marss$states[1,39]
q = coef(fit.marss)$Q
u = coef(fit.marss)$U
#next question asks for pop size in 2007 so nforeward=10
nsim=1000
nforeward = 10
#each row holds a simulation
x=matrix(NA, nsim, nforeward)
x[,1]=x0+u+rnorm(nsim,0,sqrt(q))
for(t in 2:nforeward) x[,t]=x[,t-1]+u+rnorm(nsim,0,sqrt(q))


###################################################
### code chunk number 34: hw4-prob
###################################################
#I just want the fraction of simulations that were 50,000 or above in 2007
xthresh = log(50000)
sum(x[,10]<=xthresh)/nsim


###################################################
### code chunk number 35: hw5-fit1
###################################################
mod.list=list(
  B=matrix(1), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  x0=matrix("mu"), tinitx=0)
fit.whales1 = MARSS(dat, model=mod.list)


###################################################
### code chunk number 36: hw5-fit2
###################################################
mod.list=list(
  B=matrix(1), U=matrix(0), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0),
  x0=matrix("mu"), tinitx=0)
fit.whales2 = MARSS(dat, model=mod.list)


###################################################
### code chunk number 37: hw5-fit3
###################################################
mod.list=list(
  B=matrix(1), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix(0.05),
  x0=matrix("mu"), tinitx=0)
fit.whales3 = MARSS(dat, model=mod.list)


###################################################
### code chunk number 38: hw5-aicc
###################################################
c(fit.whales1$AICc, fit.whales2$AICc, fit.whales3$AICc)
c(fit.whales1$logLik, fit.whales2$logLik, fit.whales3$logLik)


###################################################
### code chunk number 39: hw5-aicc-table
###################################################
AICc=c(fit.whales1$AICc, fit.whales2$AICc, fit.whales3$AICc)
delAIC = AICc-min(AICc)
relLik = exp(-0.5*delAIC)
aic.table=data.frame(
  AICc = AICc,
  delAICc = delAIC,
  relLik = relLik/sum(relLik)
)
rownames(aic.table) = c(
  "proc only with drift", 
  "proc only no drift", 
  "proc with drift and obs error")
round(aic.table, digits=3)


###################################################
### code chunk number 40: hw5-data
###################################################
require(forecast)
dat=log(airmiles)
n = length(dat)
training.dat = dat[1:(length(airmiles)-3)]
test.dat = dat[(length(airmiles)-2):n]


###################################################
### code chunk number 41: hw6-fit
###################################################
fit.1=Arima(training.dat, order =c(0,0,0))
fit.2=Arima(training.dat, order =c(1,0,0))
fit.3=Arima(training.dat, order =c(0,0,1))
fit.4=Arima(training.dat, order =c(1,0,1))


###################################################
### code chunk number 42: hw6-forecast
###################################################
forecast.1=forecast(fit.1, h=3)
forecast.2=forecast(fit.2, h=3)
forecast.3=forecast(fit.3, h=3)
forecast.4=forecast(fit.4, h=3)
forecast.1


###################################################
### code chunk number 43: hw-accuracy-1
###################################################
accuracy(forecast.1, test.dat)


###################################################
### code chunk number 44: hw6-table
###################################################
MASEs = c(
  accuracy(forecast.1, test.dat)["Test set","MASE"],
  accuracy(forecast.2, test.dat)["Test set","MASE"],
  accuracy(forecast.3, test.dat)["Test set","MASE"],
  accuracy(forecast.4, test.dat)["Test set","MASE"]
)
data.frame(
  name=paste("Arima",c("(0,0,0)","(1,0,0)","(0,0,1)","(1,0,1)"),sep=""), 
  MASE=MASEs
  )


###################################################
### code chunk number 45: reset
###################################################
options(prompt="> ", continue=" +", width=120)


