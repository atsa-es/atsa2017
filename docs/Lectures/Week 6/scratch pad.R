

library(MARSS)

yy <- Nile

BB <- matrix(1)
UU <- matrix(0)
QQ <- matrix("q")

ZZ <- matrix(1)
AA <- matrix(0)
RR <- matrix("r")

mod.list <- list(B=BB, U=UU, Q=QQ, Z=ZZ, A=AA, R=RR)

m1 <- MARSS(matrix(yy, nrow=1), model=mod.list)

m1up95 <- m1$states + 1.96*m1$states.se
m1lo95 <- m1$states - 1.96*m1$states.se

dev.new(height=3.3, width=5.5)

# png("Nile_DLM.png", height=3.3, width=5.5, units="in", res=300)

par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))

plot.ts(yy, type="o", pch=16, ylab="Flow of the Nile River", xlab="Year")

lines(seq(tsp(yy)[1],tsp(yy)[2]),t(m1$states), col="blue", lwd=2)
lines(seq(tsp(yy)[1],tsp(yy)[2]),t(m1lo95), col="blue", lty="dashed")
lines(seq(tsp(yy)[1],tsp(yy)[2]),t(m1up95), col="blue", lty="dashed")

# dev.off()


# trying stochstic trend/bias

BB <- matrix(c(1,0,1,1),2,2)
UU <- matrix(0,2,1)
QQ <- matrix(list(0),2,2)
diag(QQ) <- c("q.lvl","q.tnd")
diag(QQ) <- c("q","q")

ZZ <- matrix(c(1,0),1,2)
AA <- matrix(0)
RR <- matrix("r")

con.list <- list(maxit=1500)

mod.list <- list(B=BB, U=UU, Q=QQ, Z=ZZ, A=AA, R=RR)

# only need starting values for regr parameters
inits.list = list(x0=matrix(c(0, 0), nrow=2))

dat <- matrix((yy-mean(yy))/sqrt(var(yy)), nrow=1)
m2 <- MARSS(dat, model=mod.list, control=con.list, inits=inits.list)

m2up95 <- m2$states + 1.96*m2$states.se
m2lo95 <- m2$states - 1.96*m2$states.se

dev.new(height=3.3, width=5.5)

# png("Nile_DLM.png", height=3.3, width=5.5, units="in", res=300)

par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))

plot.ts(yy, type="o", pch=16, ylab="Flow of the Nile River", xlab="Year")

lines(seq(tsp(yy)[1],tsp(yy)[2]),t(m2$states), col="blue", lwd=2)
lines(seq(tsp(yy)[1],tsp(yy)[2]),t(m2lo95), col="blue", lty="dashed")
lines(seq(tsp(yy)[1],tsp(yy)[2]),t(m2up95), col="blue", lty="dashed")





