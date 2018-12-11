###################################################
### code chunk number 2: Intro to jags.Rnw:20-22
###################################################
require(coda)
require(R2jags)


###################################################
### code chunk number 3: Intro to jags.Rnw:29-33
###################################################
data(airquality)
Wind = airquality$Wind # wind speed
Temp = airquality$Temp # air temperature
N = dim(airquality)[1] # number of data points


###################################################
### code chunk number 4: Intro to jags.Rnw:49-74
###################################################
#################################################################################
# 1. START WITH AN EXAMPLE OF LINEAR REGRESSION
# no covariates, so intercept only. The parameters here are the mean 'mu' and 
# precision/variance parameter 'tau.obs'
#################################################################################
model.loc="lm_intercept.txt" # name of the txt file
jagsscript = cat("
model {  
   # priors on parameters
   mu ~ dnorm(0, 0.01); # This is normal with mean = 0, sd = 1/sqrt(0.01)
   tau.obs ~ dgamma(0.001,0.001); # This is inverse gamma
   sd.obs <- 1/sqrt(tau.obs); # sd is treated as derived parameter

   # Jags is not vectorized, so we have to loop over observations
   for(i in 1:N) {
      # for each observation, we'll assume the data to be 
      # normally distributed around a common mean
      predY[i] <- mu; 
      Y[i] ~ dnorm(predY[i], tau.obs);
      # The above 2 lines are inefficient, and could also be written as
      # Y[i] ~ dnorm(mu, tau);
   }
}  

",file=model.loc)


###################################################
### code chunk number 5: Intro to jags.Rnw:78-79 (eval = FALSE)
###################################################
## ?jags 


###################################################
### code chunk number 6: Intro to jags.Rnw:84-85 (eval = FALSE)
###################################################
## ?jags.parallel


###################################################
### code chunk number 7: Intro to jags.Rnw:90-94
###################################################
jags.data = list("Y"=Wind,"N"=N) # named list of data
jags.params=c("sd.obs","mu") # parameters in the linear regression model
mod_lm_intercept = jags(jags.data, parameters.to.save=jags.params, 
model.file=model.loc, n.chains = 3, n.burnin=5000, n.thin=1, n.iter=10000, DIC=TRUE)  


###################################################
### code chunk number 8: Intro to jags.Rnw:100-101
###################################################
mod_lm_intercept


###################################################
### code chunk number 9: Intro to jags.Rnw:106-107
###################################################
attach.jags(mod_lm_intercept)


###################################################
### code chunk number 10: Intro to jags.Rnw:112-116
###################################################
# Now we can make plots of posterior values
par(mfrow = c(2,1))
hist(mu,40,col="grey",xlab="Mean",main="")
hist(sd.obs,40,col="grey",xlab=expression(sigma[obs]),main="")


###################################################
### code chunk number 11: Intro to jags.Rnw:121-128
###################################################
createMcmcList = function(jagsmodel) {
McmcArray = as.array(jagsmodel$BUGSoutput$sims.array)
McmcList = vector("list",length=dim(McmcArray)[2])
for(i in 1:length(McmcList)) McmcList[[i]] = as.mcmc(McmcArray[,i,])
McmcList = mcmc.list(McmcList)
return(McmcList)
}


###################################################
### code chunk number 12: Intro to jags.Rnw:140-143
###################################################
myList = createMcmcList(mod_lm_intercept)
summary(myList[[1]])
plot(myList[[1]])


###################################################
### code chunk number 13: Intro to jags.Rnw:149-155 (eval = FALSE)
###################################################
## # Run the majority of the diagnostics that CODA() offers
## library(coda)
## gelmanDiags = gelman.diag(createMcmcList(mod_lm_intercept),multivariate=F)
## autocorDiags = autocorr.diag(createMcmcList(mod_lm_intercept))
## gewekeDiags = geweke.diag(createMcmcList(mod_lm_intercept))
## heidelDiags = heidel.diag(createMcmcList(mod_lm_intercept))


