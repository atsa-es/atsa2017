

tstar <- 11

x1 <-  c(rep(0,10),1,rep(0,9))

pdf("indicator_examples.pdf", height=3, width=7)

par(mai=c(0.9,0.9,0.1,0.1), omi=c(0,0,0,0))

plot.ts(x1, ylab="Indicator variable", type="o", pch=16, yaxt="n")

axis(2, at=c(0,1), las=1)

scl <- 0.07
aB <- (par()$usr[4] - par()$usr[3])*scl
arrows(x0=tstar,y0=par()$usr[3]-aB,y1=par()$usr[3], length=0.1, col="red",xpd=TRUE)
text(tstar,par()$usr[3]-aB,expression(paste(italic(t),"*")),
	 cex=1.2, col="red", pos=1, xpd=TRUE)


x2 <-  c(rep(0,10),rep(1,10))

plot.ts(x2, ylab="Indicator variable", type="o", pch=16, yaxt="n")

axis(2, at=c(0,1), las=1)

scl <- 0.07
aB <- (par()$usr[4] - par()$usr[3])*scl
arrows(x0=tstar,y0=par()$usr[3]-aB,y1=par()$usr[3], length=0.1, col="red",xpd=TRUE)
text(tstar,par()$usr[3]-aB,expression(paste(italic(t),"*")),
	 cex=1.2, col="red", pos=1, xpd=TRUE)


x3 <-  c(rep(0,10),1,0,0,-1,0,1,0,0,0)

plot.ts(x3, ylab="Indicator variable", type="o", pch=16, yaxt="n", ylim=c(-1,1))

axis(2, at=c(-1,0,1), las=1)

scl <- 0.07
aB <- (par()$usr[4] - par()$usr[3])*scl
arrows(x0=tstar,y0=par()$usr[3]-aB,y1=par()$usr[3], length=0.1, col="red",xpd=TRUE)
text(tstar,par()$usr[3]-aB,expression(paste(italic(t),"*")),
	 cex=1.2, col="red", pos=1, xpd=TRUE)
arrows(x0=14,y0=par()$usr[3]-aB,y1=par()$usr[3], length=0.1, col="red",xpd=TRUE)
text(14,par()$usr[3]-aB,"Off",
	 cex=1.2, col="red", pos=1, xpd=TRUE)
arrows(x0=16,y0=par()$usr[3]-aB,y1=par()$usr[3], length=0.1, col="red",xpd=TRUE)
text(16,par()$usr[3]-aB,"On",
	 cex=1.2, col="red", pos=1, xpd=TRUE)


x4 <- c(rep(0,10),0.94,rep(0,9))

plot.ts(x4, ylab="Expense (1000s USD)", type="o", pch=16, yaxt="n", ylim=c(0,1))

axis(2, at=c(0,1), labels=c(0,1000))

scl <- 0.07
aB <- (par()$usr[4] - par()$usr[3])*scl
arrows(x0=tstar,y0=par()$usr[3]-aB,y1=par()$usr[3], length=0.1, col="red",xpd=TRUE)
text(tstar,par()$usr[3]-aB,expression(paste(italic(t),"*")),
	 cex=1.2, col="red", pos=1, xpd=TRUE)


x5 <- c(rep(0,10),0.94,0.87,0.82,0.8,rep(0,6))

plot.ts(x5, ylab="Expense (1000s USD)", type="o", pch=16, yaxt="n", ylim=c(0,1))

axis(2, at=c(0,1), labels=c(0,1000))

scl <- 0.07
aB <- (par()$usr[4] - par()$usr[3])*scl
arrows(x0=tstar,y0=par()$usr[3]-aB,y1=par()$usr[3], length=0.1, col="red",xpd=TRUE)
text(tstar,par()$usr[3]-aB,expression(paste(italic(t),"*")),
	 cex=1.2, col="red", pos=1, xpd=TRUE)



dev.off()





