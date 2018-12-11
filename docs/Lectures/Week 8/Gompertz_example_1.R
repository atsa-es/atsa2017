#################################################################
##  Discrete Gompertz
##  Example 1.  Estimating the observation error variance is hard
#################################################################

require(MARSS)

#set.seed(123)
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

#Let's add some observation error
r=.1
a=0
y = x + rnorm(n,0,sqrt(r))+a

mod.list=list(
  U=matrix("u"),
  x0=matrix("x0"),
  B=matrix("b"),
  Q=matrix("q"),
  Z=matrix(1),
  A=matrix(0),
  R=matrix("r"),
  tinitx=1)

fit2=MARSS(y,model=mod.list,method="BFGS")

#Let's forecast our OBSERVATIONS forward
t.forward = 10

#Let's first add the the real x and observations
t0=n-1
ylims=c(x[n-t0]-5*(r+q/(1-b^2)),u/(1-b)+5*(r+q/(1-b^2)))
xlims=c(n-t0,n+t.forward)
par(mfrow=c(1,1))
plot((n-t0):n, x[(n-t0):n],xlim=xlims,ylim=ylims,type="l",ylab="y",xlab="t")
points(y)
title(paste("forecast with",n,"data points for estimation\nblue is estimate; red is true"))

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

#Let's estimate b assuming no observation error and forecast from that

mod.list=list(
  U=matrix("u"),
  x0=matrix(y[1]),
  B=matrix("b"),
  Q=matrix("q"),
  Z=matrix(1),
  A=matrix(0),
  R=matrix(0),
  tinitx=1)

#fit3=MARSS(y, model=mod.list) #for seed 100, EM algorithm has trouble
fit3=MARSS(y, model=mod.list, method="BFGS")

#Now let's forecast 1000 times using our state-space model
xend.est=fit3$states[n]
b.est=coef(fit3)$B
u.est=coef(fit3)$U
q.est=coef(fit3)$Q
r.est=coef(fit3)$R
for(i in 1:1000){
  xt=b.est*xend.est+u.est+rnorm(1,0,sqrt(q.est))
  n=40
  for(i in 2:t.forward) xt[i]=b.est*xt[i-1]+u.est+rnorm(1,0,sqrt(q.est))  
  jit=rnorm(1,0,.1)+.5
  points(n+1:t.forward+jit,xt,pch=".",col="green")
}

#fit2 is the model with process and observation error estimated
rbind(
  c(r,b,u,q,x[1]),
 coef(fit2,type="vector")
 )


#what if a were not zero?
a=-10
y2=y+a
mod.list2=mod.list
mod.list2$A=matrix("a")
fit.bad=MARSS(y2,model=mod.list2,method="BFGS")
