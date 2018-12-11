### R code from vignette source 'Key.Rnw'



###################################################
### code chunk number 2: key-set-up-data
###################################################
require(MARSS)
require(stats)
require(forecast)
# load data
#load("fitting multi ss models/plankdat.rdata")
load("plankdat.rdata")
# phytoplankton names
phytos = colnames(plankdat)[-(1:5)]
# matrix response var
dat = as.matrix(t(plankdat[,phytos]))
# length of time series
TT = dim(dat)[2]
# number of states & obs
m = n = dim(dat)[1]
# set period; data were collected monthly
period = 12


###################################################
### code chunk number 3: model-base
###################################################
  B = "diagonal and unequal" #given to you
  U = "zero" # we demeaned the data
  Q = "diagonal and unequal" #given to you

  # and for the observation eqn:
  Z = "identity" #each row of data is one species
  A = "zero" #demeaned the data
  R = "diagonal and equal" # given to you


###################################################
### code chunk number 4: c-month-as-factor
###################################################
c.fac <- diag(period)
for(i in 2:(ceiling(TT / period))) { c.fac <- cbind(c.fac, diag(period)) }
# better row names
the.months = month.abb #the months
rownames(c.fac) = the.months


###################################################
### code chunk number 5: Key.Rnw:56-60
###################################################
C = "unconstrained"; c=c.fac
D="zero"; d="zero"
model.list = list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,C=C,c=c,D=D,d=d)
q1a <- MARSS(dat, model=model.list, control=list(maxit=1500))


###################################################
### code chunk number 6: q1.fac.seas
###################################################
seas.a <- coef(q1a,type="matrix")$C
rownames(seas.a) <- phytos
colnames(seas.a) <- the.months


###################################################
### code chunk number 7: Key.Rnw:71-74
###################################################
matplot(t(seas.a),type="l",bty="n",xaxt="n", ylab="Fixed monthly", col=1:5)
axis(1,labels=the.months, at=1:12,las=2,cex.axis=0.75)
legend("topright", lty=1:5, legend=phytos, cex=0.6, col=1:5)


###################################################
### code chunk number 8: Key.Rnw:85-90
###################################################
poly.order = 3
month.cov = matrix(1,1,period)
for(i in 1:poly.order) { month.cov <- rbind(month.cov,(1:12)^i)  }
# for c, month.cov is replicated 10 times (once for each year)
c.m.poly <- matrix(month.cov, poly.order+1, TT, byrow=FALSE)


###################################################
### code chunk number 9: Key.Rnw:94-98
###################################################
C = "unconstrained"; c=c.m.poly
D="zero"; d="zero"
model.list = list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,C=C,c=c,D=D,d=d)
q1b <- MARSS(dat, model=model.list, control=list(maxit=1500))


###################################################
### code chunk number 10: Key.Rnw:102-106
###################################################
C.b = coef(q1b,type="matrix")$C
seas.b <- C.b %*% month.cov
rownames(seas.b) <- phytos
colnames(seas.b) <- the.months


###################################################
### code chunk number 11: Key.Rnw:113-116
###################################################
cos.t <- cos(2 * pi * seq(TT) / period)
sin.t <- sin(2 * pi * seq(TT) / period)
c.Four <- rbind(cos.t,sin.t)


###################################################
### code chunk number 12: Key.Rnw:120-124
###################################################
C = "unconstrained"; c=c.Four
D="zero"; d="zero"
model.list = list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,C=C,c=c,D=D,d=d)
q1c = MARSS(dat, model=model.list, control=list(maxit=1500))


###################################################
### code chunk number 13: Key.Rnw:128-132
###################################################
C.c = coef(q1c,type="matrix")$C
seas.c = C.c %*% c.Four[,1:period]
rownames(seas.c) = phytos
colnames(seas.c) = the.months


###################################################
### code chunk number 14: q1-fig-plot
###################################################
par(mfrow=c(3,1), mar=c(2,4,2,2)) 
matplot(t(seas.a),type="l",bty="n",xaxt="n", ylab="Fixed monthly", col=1:5)
axis(1,labels=the.months, at=1:12,las=2,cex.axis=0.75)

matplot(t(seas.b),type="l",bty="n",xaxt="n", ylab="Fixed monthly", col=1:5)
axis(1,labels=the.months, at=1:12,las=2,cex.axis=0.75)

matplot(t(seas.c),type="l",bty="n",xaxt="n", ylab="Fixed monthly", col=1:5)
axis(1,labels=the.months, at=1:12,las=2,cex.axis=0.75)
legend("topright", lty=1:5, legend=phytos, cex=0.6, col=1:5)


###################################################
### code chunk number 15: Key.Rnw:165-174
###################################################
covars <- rbind(Temp=plankdat[,"Temp"], TP=plankdat[,"TP"])
# z-scored covariates
prec = solve(diag(sqrt(apply(covars, 1, var))))
avg = apply(covars, 1, mean)
covars.z <-  prec %*% (covars-avg)
rownames(covars.z) <- rownames(covars)
# Ensure that mean and variance are indeed 0 and 1 (it's easy to make mistakes)
apply(covars.z,1,mean)  # Note: 1e-17 is basically 0
apply(covars.z,1,var)


###################################################
### code chunk number 16: Key.Rnw:178-182
###################################################
C = "unconstrained"
c=covars.z[1,6:TT]
for(i in 1:5) c = rbind(c,
  covars.z[1,(6-i):(TT-i)] )
rownames(c)=paste("lag",0:5,sep="")
model.list <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,C=C,c=c,D=D,d=d)
q2.Temp = MARSS(dat[,6:TT], model=model.list, control=list(maxit=1500))


###################################################
### code chunk number 17: Key.Rnw:186-190
###################################################
C = matrix(list(0),5,2); C[,2]=paste("TP",1:5,sep="")
c = covars.z[,1:(TT-1)]                         
model.list <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,C=C,c=c,D=D,d=d)
q2.TP = MARSS(dat[,2:TT], model=model.list, control=list(maxit=1500))


###################################################
### code chunk number 18: Key.Rnw:194-198
###################################################
C = "unconstrained"
c = covars.z                         
model.list <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,C=C,c=c,D=D,d=d)
q2.both = MARSS(dat, model=model.list, control=list(maxit=1500))


###################################################
### code chunk number 19: Key.Rnw:202-204
###################################################
data.frame(Model=c("Temp", "TP", "Both"),
  	   AICc=round(c(q2.Temp$AICc, q2.TP$AICc, q2.both$AICc), 1))


###################################################
### code chunk number 20: Key.Rnw:215-216
###################################################
round(coef(q2.both, type="matrix")$C)


###################################################
### code chunk number 21: Key.Rnw:224-229
###################################################
mod.names <- c("Temp", "TP", "Both", "Fixed", "Cubic", "Fourier")
AICc.all <- c(q2.Temp$AICc,q2.TP$AICc,q2.both$AICc,q1a$AICc,q1b$AICc,q1c$AIC)
delta.AICc <- round(AICc.all - min(AICc.all), 1)
# model selection results sorted from best (top) to worst (bottom)
data.frame(Model=mod.names, delta.AICc=delta.AICc)[order(delta.AICc),]


