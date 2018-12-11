### R code from vignette source 'fitting_multivariate_ss.xRnw'

###################################################
### code chunk number 1: RUNFIRST
###################################################
require(MARSS)


###################################################
### code chunk number 2: noshowlegend
###################################################
d=harborSealWA
legendnames = (unlist(dimnames(d)[2]))[2:ncol(d)]
for(i in 1:length(legendnames)) cat(paste(i,legendnames[i],"\n",sep=" "))


###################################################
### code chunk number 3: fig1
###################################################
d=harborSealWA
dat = d[,2:ncol(d)] #first col is years
x = d[,1] #first col is years
n = ncol(dat) #num time series

#set up the graphical parameters to give each data a unique line, color and width
options(warn=-99)
ltys=matrix(1,nrow=n)
cols=matrix(1:4,nrow=n)
lwds=matrix(1:2,nrow=n)
pchs=matrix(as.character(c(1:n)),nrow=n)
options(warn=0)

matplot(x,dat,xlab="",ylab="log(counts)",type="b",pch=pchs,lty=ltys,col=cols,lwd=lwds,bty="L")
title("Puget Sound Harbor Seal Surveys")


###################################################
### code chunk number 4: Cs2-showdata
###################################################
print(harborSealWA[1:8,], digits=3)


###################################################
### code chunk number 5: Cs2-readindata
###################################################
years = harborSealWA[,"Year"]
dat= harborSealWA[,!(colnames(harborSealWA) %in% c("Year", "HC"))]
dat=t(dat) #transpose to have years across columns
colnames(dat) = years
n = nrow(dat)-1


###################################################
### code chunk number 6: fit.0.model
###################################################
mod.list.0 = list(
B=matrix(1),
U=matrix("u"),
Q=matrix("q"),
Z=matrix(1,4,1),
A="scaling",
R="diagonal and unequal",
x0=matrix("mu"),
tinitx=0 )


###################################################
### code chunk number 7: fit.0.fit
###################################################
fit.0 = MARSS(dat, model=mod.list.0)


###################################################
### code chunk number 8: mss_model_resids
###################################################
par(mfrow=c(2,2))
resids=residuals(fit.0)
for(i in 1:4){
plot(resids$model.residuals[i,],ylab="model residuals", xlab="")
abline(h=0)
title(rownames(dat)[i])
}


###################################################
### code chunk number 9: mss_model_resids_plot
###################################################
par(mfrow=c(2,2))
resids=residuals(fit.0)
for(i in 1:4){
plot(resids$model.residuals[i,],ylab="model residuals", xlab="")
abline(h=0)
title(rownames(dat)[i])
}


###################################################
### code chunk number 10: fit.1.model
###################################################
mod.list.1 = list(
B="identity",
U="equal",
Q="equalvarcov",
Z="identity",
A="scaling",
R="diagonal and unequal",
x0="unequal",
tinitx=0 )


###################################################
### code chunk number 11: fit.1.fit
###################################################
fit.1 = MARSS(dat, model=mod.list.1)


###################################################
### code chunk number 12: fits.aicc
###################################################
c(fit.0$AICc, fit.1$AICc)


###################################################
### code chunk number 13: mss_model_resids_1
###################################################
par(mfrow=c(2,2))
resids=residuals(fit.1)
for(i in 1:4){
plot(resids$model.residuals[i,],ylab="model residuals", xlab="")
abline(h=0)
title(rownames(dat)[i])
}


###################################################
### code chunk number 14: fig2
###################################################
par(mfrow=c(2,2))
for(i in 1:4){
plot(years,fit.1$states[i,],ylab="log subpopulation estimate", xlab="", type="l")
lines(years,fit.1$states[i,]-1.96*fit.1$states.se[i,],type="l",lwd=1,lty=2,col="red")
lines(years,fit.1$states[i,]+1.96*fit.1$states.se[i,],type="l",lwd=1,lty=2,col="red")
title(rownames(dat)[i])
}


###################################################
### code chunk number 15: fig2_plot
###################################################
par(mfrow=c(2,2))
for(i in 1:4){
plot(years,fit.1$states[i,],ylab="log subpopulation estimate", xlab="", type="l")
lines(years,fit.1$states[i,]-1.96*fit.1$states.se[i,],type="l",lwd=1,lty=2,col="red")
lines(years,fit.1$states[i,]+1.96*fit.1$states.se[i,],type="l",lwd=1,lty=2,col="red")
title(rownames(dat)[i])
}


###################################################
### code chunk number 16: Cs01_set.up.data
###################################################
years = harborSeal[,"Year"]
#leave off Hood Canal data for now
good = !(colnames(harborSeal)%in%c("Year","HoodCanal"))
sealData = t(harborSeal[,good])


###################################################
### code chunk number 17: Cs02_fig1
###################################################
par(mfrow=c(4,3),mar=c(2,2,2,2))
for(i in 2:dim(harborSeal)[2]) {
    plot(years, harborSeal[,i], xlab="", ylab="", main=colnames(harborSeal)[i])
}


###################################################
### code chunk number 18: Z.model
###################################################
Z.model=matrix(0,11,3)
Z.model[c(1,2,9,10),1]=1  #which elements in col 1 are 1
Z.model[c(3:6,11),2]=1  #which elements in col 2 are 1
Z.model[7:8,3]=1  #which elements in col 3 are 1


###################################################
### code chunk number 19: Z.model.1
###################################################
Z1=factor(c("wa.or","wa.or",rep("ps",4),"ca","ca","wa.or","wa.or","bc")) 


###################################################
### code chunk number 20: model_list
###################################################
mod.list = list(
B = "identity",
U = "unequal",
Q = "equalvarcov",
Z = "placeholder",
A = "scaling",
R = "diagonal and equal",
x0 = "unequal",
tinitx = 0 )


###################################################
### code chunk number 21: set-up-Zs
###################################################
Z.models = list(
H1=factor(c("wa.or","wa.or",rep("ps",4),"ca","ca","wa.or","wa.or","bc")), 
H2=factor(c(rep("coast",2),rep("ps",4),rep("coast",4),"ps")), 
H3=factor(c(rep("N",6),"S","S","N","S","N")),
H4=factor(c("nc","nc","is","is","ps","ps","sc","sc","nc","sc","is")),
H5=factor(rep("pan",11)),
H6=factor(1:11) #site
)
names(Z.models)=
     c("stock","coast+PS","N+S","NC+strait+PS+SC","panmictic","site")


###################################################
### code chunk number 22: Cs05_run.the.models
###################################################
out.tab=NULL
fits=list()
for(i in 1:length(Z.models)){
     mod.list$Z = Z.models[[i]] 
     fit = MARSS(sealData, model=mod.list,
            silent=TRUE, control=list(maxit=1000))
     out=data.frame(H=names(Z.models)[i], 
            logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
            m=length(unique(Z.models[[i]])),
            num.iter=fit$numIter, converged=!fit$convergence)
     out.tab=rbind(out.tab,out)
     fits=c(fits,list(fit))
}


###################################################
### code chunk number 23: Cs06_sort.results
###################################################
min.AICc=order(out.tab$AICc)
out.tab.1=out.tab[min.AICc,]


###################################################
### code chunk number 24: Cs07_add.delta.aicc
###################################################
out.tab.1=cbind(out.tab.1,
           delta.AICc=out.tab.1$AICc-out.tab.1$AICc[1])


###################################################
### code chunk number 25: Cs08_add.delta.aicc
###################################################
out.tab.1=cbind(out.tab.1, 
           rel.like=exp(-1*out.tab.1$delta.AICc/2))


###################################################
### code chunk number 26: Cs09_aic.weight
###################################################
out.tab.1=cbind(out.tab.1,
          AIC.weight = out.tab.1$rel.like/sum(out.tab.1$rel.like))


###################################################
### code chunk number 27: Cs10_print.table
###################################################
out.tab.1$delta.AICc = round(out.tab.1$delta.AICc, digits=2)
out.tab.1$AIC.weight = round(out.tab.1$AIC.weight, digits=3)
print(out.tab.1[,c("H","delta.AICc","AIC.weight", "converged")], row.names=FALSE)



