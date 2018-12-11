#################################################################
##  Discrete Gompertz
##  Example 3.  Non-stationary Mean-reverting random walks; ML at edge
##  Finding the interior ML using sims
#################################################################

#Same code as before but different random seed and bigger q
set.seed(120)

##A mean-reverting model
K=log(2000)
u=1
b=1-u/K
x0=1
q=.2 #increased

#generate and underlying true process
x=b*x0+u+rnorm(1,0,sqrt(q))

n=40
for(i in 2:n) x[i]=b*x[i-1]+u+rnorm(1,0,sqrt(q))

#add some observation error
r=.1
y = x + rnorm(n,0,sqrt(r))

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
fit2b=MARSS(y,model=mod.list,method="BFGS")

#Look for an interior local MLE
Rlim=c(.001,seq(.01,.4,.02))
ests=data.frame(r=Rlim,LL=NA)
for(r in Rlim){
  tmp=MARSS(y,model=list(B="unconstrained",R=matrix(r),tinitx=1),silent=TRUE)
  ests$LL[ests$r==r]=tmp$logLik
  cat(r," ")
}
plot(ests$r,ests$LL,xlab="r",ylab="LL", type="l",lty=1)
r.interior=ests$r[ests$LL==max(ests$LL)]
abline(v=r.interior,col="red")
title("Log Likelihood plotted against r values")

MARSS(y,model=list(B="unconstrained",R=matrix(r.interior),tinitx=1))

readline("Continue?")

# Example using redstart data from Dennis et al 2006
# This is an estimate where the time series is not assumed to be stationary
redstart=c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2, 4,11,11,9,6) 
y=log(redstart)
plot(y, main="redstart counts")
fit.red=MARSS(y,model=list(B="unconstrained",tinitx=1),method="BFGS")

par(mfrow=c(2,1))
Rlim=seq(.01,.5,.02)
ests=data.frame(r=Rlim,LL=NA)
for(r in Rlim){
  tmp=MARSS(y,model=list(B="unconstrained",R=matrix(r),tinitx=1),method="BFGS",silent=TRUE)
  ests$LL[ests$r==r]=tmp$logLik
  cat(r," ")
}
plot(ests$r,ests$LL,xlab="r",ylab="LL", type="l",lty=1)
r.interior=ests$r[ests$LL==max(ests$LL[ests$r>.2])]
abline(v=r.interior,col="red")
title("Log Likelihood plotted against r values")

par(mfrow=c(2,1))
Qlim=seq(.01,.4,.02)
ests=data.frame(q=Qlim,LL=NA)
for(q in Qlim){
  tmp=MARSS(y,model=list(B="unconstrained",Q=matrix(q),tinitx=1),method="BFGS",silent=TRUE)
  ests$LL[ests$q==q]=tmp$logLik
  cat(q," ")
}
plot(ests$q,ests$LL,xlab="r",ylab="LL", type="l",lty=1)
q.interior=ests$q[ests$LL==max(ests$LL)]
abline(v=q.interior,col="red")
title("Log Likelihood plotted against q values")


redfit=MARSS(y,model=list(B="unconstrained",R=matrix(r.interior),tinitx=1),method="BFGS")
