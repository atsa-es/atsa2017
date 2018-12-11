### R code from vignette source 'intro-ts-funcs-content.Rnw'

###################################################
### code chunk number 1: RUNFIRST
###################################################
require(xtable)
require(forecast)
require(RCurl)
require(stringr)
tabledir="../figures/"
options(prompt=" ", continue=" ", width=60)


###################################################
### code chunk number 2: CO2data
###################################################
## get CO2 data from Mauna Loa observatory
ww1 <- "ftp://aftp.cmdl.noaa.gov/products/"
ww2 <- "trends/co2/co2_mm_mlo.txt"
CO2 <- read.table(text=getURL(paste0(ww1,ww2)))[,c(1,2,5)]
## assign better column names
colnames(CO2) <- c("year","month","ppm")


###################################################
### code chunk number 3: CO2data.load (eval = FALSE)
###################################################
## CO2 <- read.csv("CO2_data.csv")


###################################################
### code chunk number 4: CO2ts
###################################################
## create a time series (ts) object from the CO2 data
co2 <- ts(data=CO2$ppm, frequency=12,
          start=c(CO2[1,"year"],CO2[1,"month"]))


###################################################
### code chunk number 5: plotdataPar1 (eval = FALSE)
###################################################
## ## plot the ts
## plot.ts(co2, ylab=expression(paste("CO"[2]," (ppm)")))


###################################################
### code chunk number 6: plotdata1
###################################################
## set the margins & text size
par(mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## plot the ts
plot.ts(co2, ylab=expression(paste("CO"[2]," (ppm)")))


###################################################
### code chunk number 7: Temp_data
###################################################
## get N Hemisphere land & ocean temperature anomalies from NOAA
ww1 <- "https://www.ncdc.noaa.gov/cag/time-series/"
ww2 <- "global/nhem/land_ocean/p12/12/1880-2014.csv"
Temp <- read.csv(text=getURL(paste0(ww1,ww2)), skip=3)
## create ts object
tmp <- ts(data=Temp$Value, frequency=12, start=c(1880,1))


###################################################
### code chunk number 8: Temp_data_offline (eval = FALSE)
###################################################
## ## load N Hemisphere land & ocean temperature anomalies
## Temp <- read.csv("Temp_data.csv")
## ## create ts object
## tmp <- ts(data=Temp$Value, frequency=12, start=c(1880,1))


###################################################
### code chunk number 9: alignData
###################################################
## intersection (only overlapping times)
datI <- ts.intersect(co2,tmp)
## dimensions of common-time data
dim(datI)
## union (all times)
datU <- ts.union(co2,tmp)
## dimensions of all-time data
dim(datU)


###################################################
### code chunk number 10: plotdataPar2 (eval = FALSE)
###################################################
## ## plot the ts
## plot(datI, main="", yax.flip=TRUE)


###################################################
### code chunk number 11: plotdata2
###################################################
## set the margins & text size
par(mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## plot the ts
plot(datI, main="", yax.flip=TRUE)


###################################################
### code chunk number 12: makeFilter
###################################################
## weights for moving avg
fltr <- c(1/2,rep(1,times=11),1/2)/12


###################################################
### code chunk number 13: plotTrendTSa (eval = FALSE)
###################################################
## ## estimate of trend
## co2.trend <- filter(co2, filter=fltr, method="convo", sides=2)
## ## plot the trend
## plot.ts(co2.trend, ylab="Trend", cex=1)


###################################################
### code chunk number 14: plotTrendTSb
###################################################
## estimate of trend
co2.trend <- filter(co2, filter=fltr, method="convo", sides=2)
## set the margins & text size
par(mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## plot the ts
plot.ts(co2.trend, ylab="Trend", cex=1)


###################################################
### code chunk number 15: getSeason
###################################################
## seasonal effect over time
co2.1T <- co2 - co2.trend


###################################################
### code chunk number 16: plotSeasTSa (eval = FALSE)
###################################################
## ## plot the monthly seasonal effects
## plot.ts(co2.1T, ylab="Seasonal effect", xlab="Month", cex=1)


###################################################
### code chunk number 17: plotSeasTSb
###################################################
## set the margins & text size
par(mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## plot the ts
plot.ts(co2.1T, ylab="Seasonal effect plus errors", xlab="Month", cex=1)


###################################################
### code chunk number 18: getSeasonTS
###################################################
## length of ts
ll <- length(co2.1T)
## frequency (ie, 12)
ff <- frequency(co2.1T)
## number of periods (years); %/% is integer division
periods <- ll %/% ff
## index of cumulative month
index <- seq(1,ll,by=ff) - 1
## get mean by month
mm <- numeric(ff)
for(i in 1:ff) {
  mm[i] <- mean(co2.1T[index+i], na.rm=TRUE)
}
## subtract mean to make overall mean=0
mm <- mm - mean(mm)


###################################################
### code chunk number 19: plotdataPar3 (eval = FALSE)
###################################################
## ## plot the monthly seasonal effects
## plot.ts(mm, ylab="Seasonal effect", xlab="Month", cex=1)


###################################################
### code chunk number 20: plotSeasMean
###################################################
## set the margins & text size
par(mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## plot the ts
plot.ts(mm, ylab="Seasonal effect", xlab="Month", cex=1)


###################################################
### code chunk number 21: getSeasonMean
###################################################
## create ts object for season
co2.seas <- ts(rep(mm, periods+1)[seq(ll)],
               start=start(co2.1T), 
               frequency=ff)


###################################################
### code chunk number 22: getError
###################################################
## random errors over time
co2.err <- co2 - co2.trend - co2.seas


###################################################
### code chunk number 23: plotdataPar4 (eval = FALSE)
###################################################
## ## plot the obs ts, trend & seasonal effect
## plot(cbind(co2,co2.trend,co2.seas,co2.err),main="",yax.flip=TRUE)


###################################################
### code chunk number 24: plotTrSeas
###################################################
## set the margins & text size
par(mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## plot the ts
plot(cbind(co2,co2.trend,co2.seas,co2.err), main="", yax.flip=TRUE)


###################################################
### code chunk number 25: decompCO2
###################################################
## decomposition of CO2 data
co2.decomp <- decompose(co2)


###################################################
### code chunk number 26: plotDecompA (eval = FALSE)
###################################################
## ## plot the obs ts, trend & seasonal effect
## plot(co2.decomp, yax.flip=TRUE)


###################################################
### code chunk number 27: plotDecompB
###################################################
## set the margins & text size
par(mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## plot the ts
plot(co2.decomp, yax.flip=TRUE)


###################################################
### code chunk number 28: plotCO2diff2Echo (eval = FALSE)
###################################################
## ## twice-difference the CO2 data
## co2.D2 <- diff(co2, differences=2)
## ## plot the differenced data
## plot(co2.D2, ylab=expression(paste(nabla^2,"CO"[2])))


###################################################
### code chunk number 29: plotCO2diff2eval
###################################################
## twice-difference the CO2 data
co2.D2 <- diff(co2, differences=2)


###################################################
### code chunk number 30: plotCO2diff2
###################################################
## set the margins & text size
par(mar=c(4,4.5,1,1), oma=c(0,0,0,0), cex=1)
## plot the differenced data
plot(co2.D2, ylab=expression(paste(nabla^2,"CO"[2])))


###################################################
### code chunk number 31: plotCO2diff12Echo (eval = FALSE)
###################################################
## ## difference the differenced CO2 data
## co2.D2D12 <- diff(co2.D2, lag=12)
## ## plot the newly differenced data
## plot(co2.D2D12,
##      ylab=expression(paste(nabla,"(",nabla^2,"CO"[2],")")))


###################################################
### code chunk number 32: plotCO2diff12eval
###################################################
## difference the differenced CO2 data
co2.D2D12 <- diff(co2.D2, lag=12)


###################################################
### code chunk number 33: plotCO2diff12
###################################################
## set the margins & text size
par(mar=c(4,4.5,1,1), oma=c(0,0,0,0), cex=1)
## plot the newly differenced data
plot(co2.D2D12, ylab=expression(paste(nabla,"(",nabla^2,"CO"[2],")")))


###################################################
### code chunk number 34: plotACFa (eval = FALSE)
###################################################
## ## correlogram of the CO2 data
## acf(co2, lag.max=36)


###################################################
### code chunk number 35: plotACFb
###################################################
## set the margins & text size
par(mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## correlogram of the CO2 data
acf(co2, lag.max=36)


###################################################
### code chunk number 36: BetterPlotACF
###################################################
plot.acf <- function(ACFobj) {
  rr <- ACFobj$acf[-1]
  kk <- length(rr)
  nn <- ACFobj$n.used
  plot(seq(kk),rr,type="h",lwd=2,yaxs="i",xaxs="i",
       ylim=c(floor(min(rr)),1),xlim=c(0,kk+1),
       xlab="Lag",ylab="Correlation",las=1)
  abline(h=-1/nn+c(-2,2)/sqrt(nn),lty="dashed",col="blue")
  abline(h=0)
}                                                                                                            


###################################################
### code chunk number 37: betterACF (eval = FALSE)
###################################################
## ## acf of the CO2 data
## co2.acf <- acf(co2, lag.max=36)
## ## correlogram of the CO2 data
## plot.acf(co2.acf)


###################################################
### code chunk number 38: DoOurACF
###################################################
## acf of the CO2 data
co2.acf <- acf(co2, lag.max=36)


###################################################
### code chunk number 39: plotbetterACF
###################################################
## set the margins & text size
par(mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## correlogram of the CO2 data
plot.acf(co2.acf)


###################################################
### code chunk number 40: LinearACFecho (eval = FALSE)
###################################################
## ## length of ts
## nn <- 100
## ## create straight line
## tt <- seq(nn)
## ## set up plot area
## par(mfrow=c(1,2))
## ## plot line
## plot.ts(tt, ylab=expression(italic(x[t])))
## ## get ACF
## line.acf <- acf(tt, plot=FALSE)
## ## plot ACF
## plot.acf(line.acf)


###################################################
### code chunk number 41: LinearACF
###################################################
## length of ts
nn <- 100
## create straight line
tt <- seq(nn)
## get ACF
line.acf <- acf(tt)


###################################################
### code chunk number 42: plotLinearACF
###################################################
## set the margins & text size
par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## plot line
plot.ts(tt, ylab=expression(italic(x[t])))
## plot ACF
plot.acf(line.acf)


###################################################
### code chunk number 43: SineACFecho (eval = FALSE)
###################################################
## ## create sine wave
## tt <- sin(2*pi*seq(nn)/12)
## ## set up plot area
## par(mfrow=c(1,2))
## ## plot line
## plot.ts(tt, ylab=expression(italic(x[t])))
## ## get ACF
## sine.acf <- acf(tt, plot=FALSE)
## ## plot ACF
## plot.acf(sine.acf)


###################################################
### code chunk number 44: SineACF
###################################################
## create sine wave
tt <- sin(2*pi*seq(nn)/12)
## get ACF
sine.acf <- acf(tt, plot=FALSE)


###################################################
### code chunk number 45: plotSineACF
###################################################
## set the margins & text size
par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## plot line
plot.ts(tt, ylab=expression(italic(x[t])))
## plot ACF
plot.acf(sine.acf)


###################################################
### code chunk number 46: SiLineACFecho (eval = FALSE)
###################################################
## ## create sine wave with trend
## tt <- sin(2*pi*seq(nn)/12) - seq(nn)/50
## ## set up plot area
## par(mfrow=c(1,2))
## ## plot line
## plot.ts(tt, ylab=expression(italic(x[t])))
## ## get ACF
## sili.acf <- acf(tt, plot=FALSE)
## ## plot ACF
## plot.acf(sili.acf)


###################################################
### code chunk number 47: SiLiACF
###################################################
## create sine wave with trend
tt <- sin(2*pi*seq(nn)/12) - seq(nn)/50
## get ACF
sili.acf <- acf(tt, plot=FALSE)


###################################################
### code chunk number 48: plotSiLiACF
###################################################
## set the margins & text size
par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## plot line
plot.ts(tt, ylab=expression(italic(x[t])))
## plot ACF
plot.acf(sili.acf)


###################################################
### code chunk number 49: plotPACFa (eval = FALSE)
###################################################
## ## PACF of the CO2 data
## pacf(co2, lag.max=36)


###################################################
### code chunk number 50: BetterPlotPACF
###################################################
plot.pacf <- function(PACFobj) {
  rr <- PACFobj$acf
  kk <- length(rr)
  nn <- PACFobj$n.used
  plot(seq(kk),rr,type="h",lwd=2,yaxs="i",xaxs="i",
       ylim=c(floor(min(rr)),1),xlim=c(0,kk+1),
       xlab="Lag",ylab="PACF",las=1)
  abline(h=-1/nn+c(-2,2)/sqrt(nn),lty="dashed",col="blue")
  abline(h=0)
}                                                                                                            


###################################################
### code chunk number 51: plotPACFb
###################################################
## set the margins & text size
par(mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## correlogram of the CO2 data
pacf(co2, lag.max=36)


###################################################
### code chunk number 52: CO2PACFecho (eval = FALSE)
###################################################
## ## PACF of the CO2 data
## co2.pacf <- pacf(co2)
## ## correlogram of the CO2 data
## plot.acf(co2.pacf)


###################################################
### code chunk number 53: LynxSunspotCCF
###################################################
## get the matching years of sunspot data
suns <- ts.intersect(lynx,sunspot.year)[,"sunspot.year"]
## get the matching lynx data
lynx <- ts.intersect(lynx,sunspot.year)[,"lynx"]


###################################################
### code chunk number 54: plotSunsLynxEcho (eval = FALSE)
###################################################
## ## plot time series
## plot(cbind(suns,lynx), yax.flip=TRUE)


###################################################
### code chunk number 55: plotSunsLynx
###################################################
## set the margins & text size
par(mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## plot the ts
plot(cbind(suns,lynx), main="", yax.flip=TRUE)


###################################################
### code chunk number 56: plotCCFa (eval = FALSE)
###################################################
## ## CCF of sunspots and lynx
## ccf(suns, log(lynx), ylab="Cross-correlation")


###################################################
### code chunk number 57: plotCCFb
###################################################
## set the margins & text size
par(mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## CCF of sunspots and lynx
ccf(suns, lynx, ylab="Cross-correlation")


###################################################
### code chunk number 58: DWNsim
###################################################
set.seed(123)
## random normal variates
GWN <- rnorm(n=100, mean=5, sd=0.2)
## random Poisson variates
PWN <- rpois(n=50, lambda=20)


###################################################
### code chunk number 59: DWNsimPlotEcho (eval = FALSE)
###################################################
## ## set up plot region
## par(mfrow=c(1,2))
## ## plot normal variates with mean
## plot.ts(GWN)
## abline(h=5, col="blue", lty="dashed")
## ## plot Poisson variates with mean
## plot.ts(PWN)
## abline(h=20, col="blue", lty="dashed")


###################################################
### code chunk number 60: plotDWNsims
###################################################
## set the margins & text size
par(mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1, mfrow=c(1,2))
## plot normal variates with mean
plot.ts(GWN)
abline(h=5, col="blue", lty="dashed")
## plot Poisson variates with mean
plot.ts(PWN)
abline(h=20, col="blue", lty="dashed")


###################################################
### code chunk number 61: DWNacfEcho (eval = FALSE)
###################################################
## ## set up plot region
## par(mfrow=c(1,2))
## ## plot normal variates with mean
## acf(GWN, main="", lag.max=20)
## ## plot Poisson variates with mean
## acf(PWN, main="", lag.max=20)


###################################################
### code chunk number 62: plotACFdwn
###################################################
## set the margins & text size
par(mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1, mfrow=c(1,2))
## plot normal variates with mean
acf(GWN, main="", lag.max=20)
## plot Poisson variates with mean
acf(PWN, main="", lag.max=20)


###################################################
### code chunk number 63: RWsim
###################################################
## set random number seed
set.seed(123)
## length of time series
TT <- 100
## initialize {x_t} and {w_t}
xx <- ww <- rnorm(n=TT, mean=0, sd=1)
## compute values 2 thru TT
for(t in 2:TT) { xx[t] <- xx[t-1] + ww[t] }


###################################################
### code chunk number 64: plotRWecho (eval = FALSE)
###################################################
## ## setup plot area
## par(mfrow=c(1,2))
## ## plot line
## plot.ts(xx, ylab=expression(italic(x[t])))
## ## plot ACF
## plot.acf(acf(xx, plot=FALSE))


###################################################
### code chunk number 65: calcRWACF
###################################################
xx.acf <- acf(xx, plot=FALSE)


###################################################
### code chunk number 66: plotRW
###################################################
## setup plot area
par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## plot line
plot.ts(xx, ylab=expression(italic(x[t])))
## plot ACF
plot.acf(xx.acf)


###################################################
### code chunk number 67: RWsimAlt
###################################################
## simulate RW
x2 <- cumsum(ww)


###################################################
### code chunk number 68: plotRWsimEcho (eval = FALSE)
###################################################
## ## setup plot area
## par(mfrow=c(1,2))
## ## plot 1st RW
## plot.ts(xx, ylab=expression(italic(x[t])))
## ## plot 2nd RW
## plot.ts(x2, ylab=expression(italic(x[t])))


###################################################
### code chunk number 69: plotRWalt
###################################################
## setup plot area
par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(0,0,0,0), cex=1)
## plot 1st RW
plot.ts(xx, ylab=expression(italic(x[t])))
## plot 2nd RW
plot.ts(x2, ylab=expression(italic(x[t])))


###################################################
### code chunk number 70: simAR1
###################################################
set.seed(456)
## list description for AR(1) model with small coef
AR.sm <- list(order=c(1,0,0), ar=0.1, sd=0.1)
## list description for AR(1) model with large coef
AR.lg <- list(order=c(1,0,0), ar=0.9, sd=0.1)
## simulate AR(1)
AR1.sm <- arima.sim(n=50, model=AR.sm)
AR1.lg <- arima.sim(n=50, model=AR.lg)


###################################################
### code chunk number 71: plotAR1sims (eval = FALSE)
###################################################
## ## setup plot region
## par(mfrow=c(1,2))
## ## get y-limits for common plots
## ylm <- c(min(AR1.sm,AR1.lg), max(AR1.sm,AR1.lg))
## ## plot the ts
## plot.ts(AR1.sm, ylim=ylm,
##         ylab=expression(italic(x)[italic(t)]),
##         main=expression(paste(phi," = 0.1")))
## plot.ts(AR1.lg, ylim=ylm,
##         ylab=expression(italic(x)[italic(t)]),
##         main=expression(paste(phi," = 0.9")))


###################################################
### code chunk number 72: getPlotLims
###################################################
## get y-limits for common plots
ylm <- c(min(AR1.sm,AR1.lg), max(AR1.sm,AR1.lg))


###################################################
### code chunk number 73: plotAR1contrast
###################################################
## set the margins & text size
par(mfrow=c(1,2), mar=c(4,4,1.5,1), oma=c(0,0,0,0), cex=1)
## plot the ts
plot.ts(AR1.sm, ylim=ylm,
        ylab=expression(italic(x)[italic(t)]),
        main=expression(paste(phi," = 0.1")))
plot.ts(AR1.lg, ylim=ylm,
        ylab=expression(italic(x)[italic(t)]),
        main=expression(paste(phi," = 0.9")))


###################################################
### code chunk number 74: simAR1opps
###################################################
set.seed(123)
## list description for AR(1) model with small coef
AR.pos <- list(order=c(1,0,0), ar=0.5, sd=0.1)
## list description for AR(1) model with large coef
AR.neg <- list(order=c(1,0,0), ar=-0.5, sd=0.1)
## simulate AR(1)
AR1.pos <- arima.sim(n=50, model=AR.pos)
AR1.neg <- arima.sim(n=50, model=AR.neg)


###################################################
### code chunk number 75: plotAR1oppsEcho (eval = FALSE)
###################################################
## ## setup plot region
## par(mfrow=c(1,2))
## ## get y-limits for common plots
## ylm <- c(min(AR1.pos,AR1.neg), max(AR1.pos,AR1.neg))
## ## plot the ts
## plot.ts(AR1.pos, ylim=ylm,
##         ylab=expression(italic(x)[italic(t)]),
##         main=expression(paste(phi[1]," = 0.5")))
## plot.ts(AR1.neg,
##         ylab=expression(italic(x)[italic(t)]),
##         main=expression(paste(phi[1]," = -0.5")))


###################################################
### code chunk number 76: getPlotLimsOpps
###################################################
## get y-limits for common plots
ylm <- c(min(AR1.pos,AR1.neg), max(AR1.pos,AR1.neg))


###################################################
### code chunk number 77: plotAR1opps
###################################################
## set the margins & text size
par(mfrow=c(1,2), mar=c(4,4,1.5,1), oma=c(0,0,0,0), cex=1)
## plot the ts
plot.ts(AR1.pos, ylim=ylm,
        ylab=expression(italic(x)[italic(t)]),
        main=expression(paste(phi[1]," = 0.5")))
plot.ts(AR1.neg, ylim=ylm,
        ylab=expression(italic(x)[italic(t)]),
        main=expression(paste(phi[1]," = -0.5")))


###################################################
### code chunk number 78: ARpFail (eval = FALSE)
###################################################
## arima.sim(n=100, model=list(order(2,0,0), ar=c(0.5,0.5)))


###################################################
### code chunk number 79: ARpSims
###################################################
set.seed(123)
## the 4 AR coefficients
ARp <- c(0.7, 0.2, -0.1, -0.3)
## empty list for storing models
AR.mods <- list()
## loop over orders of p
for(p in 1:4) {
  ## assume SD=1, so not specified
  AR.mods[[p]] <- arima.sim(n=10000, list(ar=ARp[1:p]))
}


###################################################
### code chunk number 80: plotARpCompsEcho (eval = FALSE)
###################################################
## ## set up plot region
## par(mfrow=c(4,3))
## ## loop over orders of p
## for(p in 1:4) {
##   plot.ts(AR.mods[[p]][1:50],
##           ylab=paste("AR(",p,")",sep=""))
##   acf(AR.mods[[p]], lag.max=12)
##   pacf(AR.mods[[p]], lag.max=12, ylab="PACF")
## }


###################################################
### code chunk number 81: plotARpComps
###################################################
## set the margins & text size
par(mfrow=c(4,3), mar=c(4,4,0.5,0.5), oma=c(0,0,0,0), cex=1)
## loop over orders of p
for(p in 1:4) {
  plot.ts(AR.mods[[p]][1:50],ylab=paste("AR(",p,")",sep=""))
  acf(AR.mods[[p]], lag.max=12)
  pacf(AR.mods[[p]], lag.max=12, ylab="PACF")
}


###################################################
### code chunk number 82: simMA1opps
###################################################
set.seed(123)
## list description for MA(1) model with small coef
MA.sm <- list(order=c(0,0,1), ma=0.2, sd=0.1)
## list description for MA(1) model with large coef
MA.lg <- list(order=c(0,0,1), ma=0.8, sd=0.1)
## list description for MA(1) model with large coef
MA.neg <- list(order=c(0,0,1), ma=-0.5, sd=0.1)
## simulate MA(1)
MA1.sm <- arima.sim(n=50, model=MA.sm)
MA1.lg <- arima.sim(n=50, model=MA.lg)
MA1.neg <- arima.sim(n=50, model=MA.neg)


###################################################
### code chunk number 83: plotMA1oppsEcho (eval = FALSE)
###################################################
## ## setup plot region
## par(mfrow=c(1,3))
## ## plot the ts
## plot.ts(MA1.sm,
##         ylab=expression(italic(x)[italic(t)]),
##         main=expression(paste(theta," = 0.2")))
## plot.ts(MA1.lg,
##         ylab=expression(italic(x)[italic(t)]),
##         main=expression(paste(theta," = 0.8")))
## plot.ts(MA1.neg,
##         ylab=expression(italic(x)[italic(t)]),
##         main=expression(paste(theta," = -0.5")))


###################################################
### code chunk number 84: plotMA1opps
###################################################
## set the margins & text size
par(mfrow=c(1,3), mar=c(4,4,1.5,0.5), oma=c(0,0,0,0), cex=1)
## plot the ts
plot.ts(MA1.sm,
        ylab=expression(italic(x)[italic(t)]),
        main=expression(paste(theta," = 0.2")))
plot.ts(MA1.lg,
        ylab=expression(italic(x)[italic(t)]),
        main=expression(paste(theta," = 0.8")))
plot.ts(MA1.neg,
        ylab=expression(italic(x)[italic(t)]),
        main=expression(paste(theta," = -0.5")))


###################################################
### code chunk number 85: MAqSims
###################################################
set.seed(123)
## the 4 MA coefficients
MAq <- c(0.7, 0.2, -0.1, -0.3)
## empty list for storing models
MA.mods <- list()
## loop over orders of q
for(q in 1:4) {
  ## assume SD=1, so not specified
  MA.mods[[q]] <- arima.sim(n=1000, list(ma=MAq[1:q]))
}


###################################################
### code chunk number 86: plotMApCompsEcho (eval = FALSE)
###################################################
## ## set up plot region
## par(mfrow=c(4,3))
## ## loop over orders of q
## for(q in 1:4) {
##   plot.ts(MA.mods[[q]][1:50],
##           ylab=paste("MA(",q,")",sep=""))
##   acf(MA.mods[[q]], lag.max=12)
##   pacf(MA.mods[[q]], lag.max=12, ylab="PACF")
## }


###################################################
### code chunk number 87: plotMApComps
###################################################
## set the margins & text size
par(mfrow=c(4,3), mar=c(4,4,0.5,0.5), oma=c(0,0,0,0), cex=1)
## loop over orders of q
for(q in 1:4) {
  plot.ts(MA.mods[[q]][1:50],ylab=paste("MA(",q,")",sep=""))
  acf(MA.mods[[q]], lag.max=12)
  pacf(MA.mods[[q]], lag.max=12, ylab="PACF")
}


###################################################
### code chunk number 88: ARMAest
###################################################
set.seed(123)
## ARMA(2,2) description for arim.sim()
ARMA22 <- list(order=c(2,0,2), ar=c(-0.7,0.2), ma=c(0.7,0.2))
## mean of process
mu <- 5
## simulated process (+ mean)
ARMA.sim <- arima.sim(n=10000, model=ARMA22) + mu
## estimate parameters
arima(x=ARMA.sim, order=c(2,0,2))


###################################################
### code chunk number 89: ARMAsearch1
###################################################
## empty list to store model fits
ARMA.res <- list()
## set counter
cc <- 1
## loop over AR
for(p in 0:3) {
  ## loop over MA
  for(q in 0:3) {
    ARMA.res[[cc]] <- arima(x=ARMA.sim,order=c(p,0,q))
    cc <- cc + 1
  }
}
## get AIC values for model evaluation
ARMA.AIC <- sapply(ARMA.res,function(x) x$aic)
## model with lowest AIC is the best
ARMA.res[[which(ARMA.AIC==min(ARMA.AIC))]]


###################################################
### code chunk number 90: autoARIMA
###################################################
## (install if necessary) & load forecast pkg
if(!require("forecast")) {
    install.packages("forecast")
    library("forecast")
}
## find best ARMA(p,q) model
auto.arima(ARMA.sim, start.p=0, max.p=3, start.q=0, max.q=3)


###################################################
### code chunk number 91: HW1_pre (eval = FALSE)
###################################################
## ## get phytoplankton data
## pp <- "http://faculty.washington.edu/scheuerl/phytoDat.txt"
## pDat <- read.table(pp)


###################################################
### code chunk number 92: HW1_1 (eval = FALSE)
###################################################
## ## what day of 2014 is Dec 1st?
## dBegin <- as.Date("2014-12-01")
## dayOfYear <- (dBegin - as.Date("2014-01-01") + 1)


###################################################
### code chunk number 93: reset
###################################################
options(prompt="> ", continue=" +", width=120)


