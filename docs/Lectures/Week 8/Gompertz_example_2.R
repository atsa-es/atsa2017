#################################################################
##  Discrete Gompertz
##  Example 2.  Example 1 with replication
#################################################################

require(MARSS)

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
#plot data
par(mfrow=c(3,1))
plot(1:n,exp(x),xlim=c(0,n),type="l",lwd=2)
acf(x,main="")
pacf(x,main="")

#Let's add some observation error
r=.1
y1 = x + rnorm(n,0,sqrt(r))
y2 = x + rnorm(n,0,sqrt(r))
y3 = x + rnorm(n,0,sqrt(r))
y=rbind(y1,y2,y3)

mod.list=list(
  U=matrix("u"),
  x0=matrix("x0"),
  B=matrix("b"),
  Q=matrix("q"),
  Z=matrix(1,3,1),
  A="zero",
  R="diagonal and equal",
  tinitx=1)

fit2=MARSS(y,model=mod.list, method="BFGS")

#Let's try simulating forward and comparing

#Let's forecast our OBSERVATIONS forward
t.forward = 10

#Let's first add the the real x and observations
t0=10
ylims=c(x[n-t0]-5*(r+q/(1-b^2)),u/(1-b)+5*(r+q/(1-b^2)))
xlims=c(n-t0,n+t.forward)
par(mfrow=c(1,1))
plot((n-t0):n, x[(n-t0):n],xlim=xlims,ylim=ylims,type="l",ylab="y",xlab="t")
points(y[1,],pch=1)
points(y[2,],pch=2)
points(y[3,],pch=3)
title(paste("forecast with",n,"data points for estimation\nred is true; blue is estimate"))

#Now let's forecast 1000 times using our state-space model
xend.est=fit2$states[n]
b.est=coef(fit2)$B
u.est=coef(fit2)$U
q.est=coef(fit2)$Q
r.est=coef(fit2)$R
for(i in 1:1000){
  xt=b.est*xend.est+u.est+rnorm(1,0,sqrt(q.est))
  n=40
  for(i in 2:t.forward) xt[i]=b.est*xt[i-1]+u.est+rnorm(1,0,sqrt(q.est))  
  jit=rnorm(1,0,.1)-0
  points(n+1:t.forward+jit,xt,pch=".",col="blue")
}

#Now let's forecast 1000 times using the truth
for(i in 1:1000){
  xt=b*x[n]+u+rnorm(1,0,sqrt(q))
  n=40
  for(i in 2:t.forward) xt[i]=b*xt[i-1]+u+rnorm(1,0,sqrt(q))  
  jit=rnorm(1,0,.1)+.25
  points(n+1:t.forward+jit,xt,pch=".",col="red")
}

