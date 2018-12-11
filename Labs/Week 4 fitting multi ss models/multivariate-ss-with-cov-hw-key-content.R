### R code from vignette source 'multivariate-ss-with-cov-hw-key-content.Rnw'

###################################################
### code chunk number 1: RUNFIRST
###################################################
options(prompt=" ", continue=" ", width=60)


###################################################
### code chunk number 2: hw4-key-set-up-data
###################################################
require(MARSS)
require(forecast)
phytos = c("Cryptomonas", "Diatoms", "Greens",
            "Unicells", "Other.algae")
yrs = lakeWAplanktonTrans[,"Year"]%in%1985:1994
dat = t(lakeWAplanktonTrans[yrs,phytos])
#z-score the data
avg = apply(dat, 1, mean, na.rm=TRUE)
sd = sqrt(apply(dat, 1, var, na.rm=TRUE))
dat = (dat-avg)/sd
rownames(dat)=phytos
#z-score the covariates
covars = rbind(Temp=lakeWAplanktonTrans[yrs,"Temp"],
               TP=lakeWAplanktonTrans[yrs,"TP"])
avg = apply(covars, 1, mean)
sd = sqrt(apply(covars, 1, var, na.rm=TRUE))
covars =  (covars-avg)/sd
rownames(covars) = c("Temp","TP")
#
#always check that the mean and variance are 1 after z-scoring
apply(dat,1,mean,na.rm=TRUE) #this should be 0
apply(dat,1,var,na.rm=TRUE) #this should be 1


###################################################
### code chunk number 3: hw4-set-up-constants
###################################################
TT = dim(dat)[2] # length of time series
m = n = dim(dat)[1] # number of states & obs
period = 12 # data were collected monthly


###################################################
### code chunk number 4: model-base
###################################################
common=list(
  B = "diagonal and unequal",
  U = "zero",
  Q = "diagonal and unequal",
  Z = "identity",
  A = "zero",
  R = "diagonal and equal",
  tinitx = 0
)
ctl = list(maxit=500) #in case we want to compare


###################################################
### code chunk number 5: q1-c-month-as-factor
###################################################
c.fac = diag(period)
for(i in 2:(ceiling(TT / period))) { c.fac = cbind(c.fac, diag(period)) }
the.months = month.abb
rownames(c.fac) = the.months


###################################################
### code chunk number 6: q1-set-up-model
###################################################
C = "unconstrained"; c=c.fac
D="zero"; d="zero"
model.list = c(common, list(C=C,c=c,D=D,d=d))
q1a = MARSS(dat, model=model.list, control=ctl)


###################################################
### code chunk number 7: q1-fac-seas
###################################################
seas.a = coef(q1a,type="matrix")$C
rownames(seas.a) = phytos
colnames(seas.a) = the.months


###################################################
### code chunk number 8: q1-plot1
###################################################
matplot(t(seas.a),type="l",bty="n",xaxt="n", ylab="Fixed monthly", col=1:5)
axis(1,labels=the.months, at=1:12,las=2,cex.axis=0.75)
legend("topright", lty=1:5, legend=phytos, cex=0.6, col=1:5)


###################################################
### code chunk number 9: q1-seas-poly
###################################################
poly.order = 3
month.cov = matrix(1,1,period)
for(i in 1:poly.order) { month.cov = rbind(month.cov,(1:12)^i)  }
# for c, month.cov is replicated 10 times (once for each year)
c.m.poly = matrix(month.cov, poly.order+1, TT, byrow=FALSE)


###################################################
### code chunk number 10: q1-seas-poly-fit
###################################################
C = "unconstrained"; c=c.m.poly
D="zero"; d="zero"
model.list = c(common, list(C=C,c=c,D=D,d=d))
q1b = MARSS(dat, model=model.list, control=ctl)


###################################################
### code chunk number 11: q1-seas-poly-get-coef
###################################################
C.b = coef(q1b,type="matrix")$C
seas.b = C.b %*% month.cov
rownames(seas.b) = phytos
colnames(seas.b) = the.months


###################################################
### code chunk number 12: q1-seas-fourier
###################################################
cos.t = cos(2 * pi * seq(TT) / period)
sin.t = sin(2 * pi * seq(TT) / period)
c.Four = rbind(cos.t,sin.t)


###################################################
### code chunk number 13: q1-seas-fourier-fit
###################################################
C = "unconstrained"; c=c.Four
D="zero"; d="zero"
model.list = c(common, list(C=C,c=c,D=D,d=d))
q1c = MARSS(dat, model=model.list, control=ctl)


###################################################
### code chunk number 14: q1-seas-fourier-coef
###################################################
C.c = coef(q1c,type="matrix")$C
seas.c = C.c %*% c.Four[,1:period]
rownames(seas.c) = phytos
colnames(seas.c) = the.months


###################################################
### code chunk number 15: q1-plot2
###################################################
par(mfrow=c(3,1), mar=c(2,4,2,2))
matplot(t(seas.a),type="l",bty="n",xaxt="n", ylab="Fixed monthly", col=1:5, lwd=2)
axis(1,labels=the.months, at=1:12,las=2,cex.axis=0.75)

matplot(t(seas.b),type="l",bty="n",xaxt="n", ylab="Fixed monthly", col=1:5, lwd=2)
axis(1,labels=the.months, at=1:12,las=2,cex.axis=0.75)

matplot(t(seas.c),type="l",bty="n",xaxt="n", ylab="Fixed monthly", col=1:5,
lwd=2)
axis(1,labels=the.months, at=1:12,las=2,cex.axis=0.75)
legend("topright", lty=1:5, legend=phytos, cex=0.6, col=1:5)


###################################################
### code chunk number 16: q2-option-1
###################################################
C = "unconstrained"
c = covars
D = d = "zero"
model.list = c(common, list(C=C,c=c,D=D,d=d))
q2.both = MARSS(dat, model=model.list, control=ctl)


###################################################
### code chunk number 17: hw2-fig1
###################################################
plot(c(rep(1,5),rep(2,5)),coef(q2.both)$C,xlab="",ylab="effect",xaxt="n",xlim=c(0,3))
axis(1, at=c(1,2),labels=c("Temp","TP"))
abline(h=0)
title("effect (C) estimates)")


###################################################
### code chunk number 18: q2-set-up
###################################################
c = covars
D = d = "zero"


###################################################
### code chunk number 19: q2-CTemp-set-up
###################################################
C = matrix(list(0),5,2); C[,1]=paste("Temp",1:5,sep="")
C


###################################################
### code chunk number 20: q2-temp-fit
###################################################
model.list = c(common, list(C=C,c=c,D=D,d=d))
q2.Temp = MARSS(dat, model=model.list, control=ctl)


###################################################
### code chunk number 21: q2-CTP-set-up
###################################################
C = matrix(list(0),5,2); C[,2]=paste("TP",1:5,sep="")
C


###################################################
### code chunk number 22: q2-TP-fit
###################################################
model.list = c(common, list(C=C,c=c,D=D,d=d))
q2.TP = MARSS(dat, model=model.list, control=ctl)


###################################################
### code chunk number 23: q2-both-fit
###################################################
C = "unconstrained"
model.list = c(common, list(C=C,c=c,D=D,d=d))
q2.both = MARSS(dat, model=model.list, control=ctl)


###################################################
### code chunk number 24: q2-dataframe
###################################################
data.frame(Model=c("Temp", "TP", "Both"),
  	   AICc=round(c(q2.Temp$AICc, q2.TP$AICc, q2.both$AICc), 1))


###################################################
### code chunk number 25: q2-model-compare
###################################################
mod.names = c("Temp", "TP", "Both", "Fixed", "Cubic", "Fourier")
AICc.all = c(q2.Temp$AICc,q2.TP$AICc,q2.both$AICc,q1a$AICc,q1b$AICc,q1c$AIC)
delta.AICc = round(AICc.all - min(AICc.all), 1)
# model selection results sorted from best (top) to worst (bottom)
data.frame(Model=mod.names, delta.AICc=delta.AICc)[order(delta.AICc),]


###################################################
### code chunk number 26: hw4-q2-fig1
###################################################
par(mfrow=c(5,2), mai=c(0.6,0.6,0.2,0.2), omi=c(0,0,0,0))
for(i in 1:5) {
  plot.ts(residuals(q2.Temp)$model.residuals[i,], ylab="Residual", main=phytos[i])
	abline(h=0, lty="dashed")
	acf(residuals(q2.Temp)$model.residuals[i,])
	}


###################################################
### code chunk number 27: hw4_q2_fig2
###################################################
par(mfrow=c(5,2), mai=c(0.6,0.6,0.2,0.2), omi=c(0,0,0,0))
for(i in 1:5) {
	plot.ts(residuals(q2.TP)$model.residuals[i,], ylab="Residual", main=phytos[i])
	abline(h=0, lty="dashed")
	acf(residuals(q2.TP)$model.residuals[i,])
	}


###################################################
### code chunk number 28: hw4-q3-set-up
###################################################
D = "unconstrained"
d = covars["Temp",,drop=FALSE]
C = c = "zero"


###################################################
### code chunk number 29: multivariate-ss-with-cov-hw-key-content.Rnw:303-305
###################################################
model.list = c(common, list(C=C,c=c,D=D,d=d))
q3.Tobs = MARSS(dat, model=model.list, control=ctl)


###################################################
### code chunk number 30: multivariate-ss-with-cov-hw-key-content.Rnw:308-310
###################################################
data.frame(Model=c("proc", "obs"),
  	   AICc=round(c(q2.Temp$AICc, q3.Tobs$AICc), 1))


###################################################
### code chunk number 31: multivariate-ss-with-cov-hw-key-content.Rnw:348-354
###################################################
# 1st: the Temp model
C = "equal"
c = covars["Temp",,drop=FALSE]
D = d = "zero"
model.list = c(common, list(C=C,c=c,D=D,d=d))
q4.Temp = MARSS(dat, model=model.list, control=ctl)


###################################################
### code chunk number 32: multivariate-ss-with-cov-hw-key-content.Rnw:357-359
###################################################
data.frame(Effect=c("taxon-specific", "equal"),
  	   AICc=round(c(q2.Temp$AICc, q4.Temp$AICc), 1))


###################################################
### code chunk number 33: multivariate-ss-with-cov-hw-key-content.Rnw:365-371
###################################################
# 2nd: the TP model
C = "equal"
c = covars["TP",,drop=FALSE]
D = d = "zero"
model.list = c(common, list(C=C,c=c,D=D,d=d))
q4.TP = MARSS(dat, model=model.list, control=ctl)


###################################################
### code chunk number 34: multivariate-ss-with-cov-hw-key-content.Rnw:374-376
###################################################
data.frame(Effect=c("taxon-specific", "equal"),
		   AICc=round(c(q2.TP$AICc, q4.TP$AICc), 1))


###################################################
### code chunk number 35: multivariate-ss-with-cov-hw-key-content.Rnw:387-393
###################################################
# 1st: unconstrained
C = "unconstrained"
c = covars
D = d = "zero"
model.list = c(common, list(C=C,c=c,D=D,d=d))
q4.unc = MARSS(dat, model=model.list, control=ctl)


###################################################
### code chunk number 36: multivariate-ss-with-cov-hw-key-content.Rnw:396-403
###################################################
# 2nd: the TP unconstrained, Temp constrained
C = matrix("Temp",5,2)
C[,2] = paste("TP",1:5)
c = covars
D = d = "zero"
model.list = c(common, list(C=C,c=c,D=D,d=d))
q4.TPunc = MARSS(dat, model=model.list, control=ctl)


###################################################
### code chunk number 37: multivariate-ss-with-cov-hw-key-content.Rnw:406-413
###################################################
# 3rd: the TP constrained, Temp unconstrained
C = matrix("TP",5,2)
C[,1] = paste("Temp",1:5)
c = covars
D = d = "zero"
model.list = c(common, list(C=C,c=c,D=D,d=d))
q4.Tempunc = MARSS(dat, model=model.list, control=ctl)


###################################################
### code chunk number 38: multivariate-ss-with-cov-hw-key-content.Rnw:416-423
###################################################
# 4th: both constrained
C = matrix("TP",5,2)
C[,1] = "Temp"
c = covars
D = d = "zero"
model.list = c(common, list(C=C,c=c,D=D,d=d))
q4.cons = MARSS(dat, model=model.list, control=ctl)


###################################################
### code chunk number 39: multivariate-ss-with-cov-hw-key-content.Rnw:425-428
###################################################
data.frame(Effect=c("unconstrained", "TP cons", "Temp cons", "both cons"),
		   AICc=round(c(q4.unc$AICc, q4.Tempunc$AICc,
		                q4.TPunc$AICc, q4.cons$AICc), 1))


###################################################
### code chunk number 40: q5-set-up
###################################################
seaslm.a = matrix(NA, 5,12, dimnames=list(phytos, month.abb))
seaslm.b = matrix(NA, 5,12, dimnames=list(phytos, month.abb))
diffdat = t(diff(t(cbind(NA,dat))))


###################################################
### code chunk number 41: q5-c-month-as-factor
###################################################
c.fac = rep(month.abb,10)
for(taxon in 1:5){
  fit=lm(dat[taxon,]~-1+c.fac)
  seaslm.a[taxon,] = coef(fit)
  fit=lm(diffdat[taxon,]~-1+c.fac)
  seaslm.b[taxon,] = coef(fit)
}


###################################################
### code chunk number 42: q5-plot1
###################################################
matplot(t(seaslm.a),type="l",bty="n",xaxt="n", ylab="Fixed monthly", col=1:5)
axis(1,labels=month.abb, at=1:12,las=2,cex.axis=0.75)
legend("topright", lty=1:5, legend=phytos, cex=0.6, col=1:5)


###################################################
### code chunk number 43: q5-plot1
###################################################
par(mfrow=c(3,1))
matplot(t(seas.a),type="l",bty="n",xaxt="n", ylab="Fixed monthly", col=1:5)
axis(1,labels=month.abb, at=1:12,las=2,cex.axis=0.75)
legend("topright", lty=1:5, legend=phytos, cex=0.6, col=1:5)
title("MARSS")

matplot(t(seaslm.a),type="l",bty="n",xaxt="n", ylab="Fixed monthly", col=1:5)
axis(1,labels=month.abb, at=1:12,las=2,cex.axis=0.75)
legend("topright", lty=1:5, legend=phytos, cex=0.6, col=1:5)
title("lm on abundance")

matplot(t(seaslm.b),type="l",bty="n",xaxt="n", ylab="Fixed monthly", col=1:5)
axis(1,labels=month.abb, at=1:12,las=2,cex.axis=0.75)
legend("topright", lty=1:5, legend=phytos, cex=0.6, col=1:5)
title("lm on month to month change")


###################################################
### code chunk number 44: q5-set-up2
###################################################
Clm = matrix(NA, 5,2, dimnames=list(phytos, c("Temp","TP")))
diffdat = t(diff(t(cbind(NA,dat))))


###################################################
### code chunk number 45: q5-c-temp-TP
###################################################
for(taxon in 1:5){
  fit=lm(dat[taxon,]~-1+covars[1,]+covars[2,])
  Clm[taxon,] = coef(fit)
}


###################################################
### code chunk number 46: hw2-fig1
###################################################
plot(c(rep(1,5),rep(2,5)),coef(q2.both)$C,xlab="",ylab="effect",xaxt="n",xlim=c(0,3),ylim=c(-1,1))
points(c(rep(1,5),rep(2,5))+.25,Clm,xlab="",col="red")
axis(1, at=c(1,2),labels=c("Temp","TP"))
abline(h=0)
title("effect (C) estimates)")
legend("topright",c("state-space","lm"),pch=1,col=c("black","red"))


###################################################
### code chunk number 47: q5-c-temp-TP
###################################################
AIClm = matrix(NA, 5,2, dimnames=list(phytos, c("Temp","TP")))
for(taxon in 1:5){
  fit=lm(dat[taxon,]~-1+covars[1,])
  AIClm[taxon,1] = AIC(fit)
  fit=lm(dat[taxon,]~-1+covars[2,])
  AIClm[taxon,2] = AIC(fit)
}
AIClm


###################################################
### code chunk number 48: q5-c-temp-tp-taxon effect
###################################################
taxon = rep(phytos,each=120)
longdat=as.vector(t(diffdat))
longtemp=rep(covars[1,],5)
longtp=rep(covars[2,],5)
fit1=lm(longdat~-1+taxon:longtemp+taxon:longtp)
fit2=lm(longdat~-1+longtemp + taxon:longtp)
fit3=lm(longdat~-1+taxon:longtemp + longtp)
fit4=lm(longdat~-1+longtemp+longtp)


###################################################
### code chunk number 49: multivariate-ss-with-cov-hw-key-content.Rnw:550-557
###################################################
data.frame(
  Effect=c("unconstrained", "TP cons", "Temp cons", "both cons"),
	AIC.lm=round(c(AIC(fit1), AIC(fit3),
		                AIC(fit2), AIC(fit4)), 1),
	AICc.ss=round(c(q4.unc$AICc, q4.Tempunc$AICc,
		                q4.TPunc$AICc, q4.cons$AICc), 1)
)


###################################################
### code chunk number 50: reset
###################################################
options(prompt="> ", continue=" +", width=120)


