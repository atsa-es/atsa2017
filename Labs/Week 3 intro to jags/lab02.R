
library(R2jags)
library(coda)

######################################################################################
# WE'LL USE SOME OF THE BUILT IN DATA IN THE DATASETS PACKAGE
# specifically 'airquality'
######################################################################################
data(airquality)
Wind = airquality$Wind
Temp = airquality$Temp
N = dim(airquality)[1]

######################################################################################
# 1. START WITH AN EXAMPLE OF LINEAR REGRESSION
# no covariates, so intercept only. The parameters here are the mean 'mu' and 
# precision/variance parameter 'tau.obs'
######################################################################################
jagsscript = cat("
model {	
   # priors on parameters
   mu ~ dnorm(0, 0.01); # This is normal with mean = 0, sd = 1/sqrt(0.01)
   tau.obs ~ dgamma(0.001,0.001); # This is inverse gamma
   sd.obs <- 1/sqrt(tau.obs); # sd is treated as derived parameter

   # Jags is not vectorized, so we have to loop over observations
   for(i in 1:N) {
      # for each observation, we'll assume the data to be normally distributed around
      # a common mean
      predY[i] <- mu; 
      Y[i] ~ dnorm(predY[i], tau.obs);
      # The above 2 lines are inefficient, and could also be written as Y[i] ~ dnorm(mu, tau);
   }
}  

",file="lm_intercept.txt")

jags.data = list("Y"=Wind,"N"=N)
jags.params=c("sd.obs","predY","mu")
model.loc=("lm_intercept.txt")
mod_lm_intercept = jags(jags.data, parameters.to.save=jags.params, model.file=model.loc, 
n.chains = 3, n.burnin=5000, n.thin=1, n.iter=10000, DIC=TRUE)  

######################################################################################
# LOOK AT SOME DIAGNOSTICS FOR THE FITTED MODEL
# 
######################################################################################

# We can attach the model to the R workspace, meaning that we have direct access to the 
# parameters we chose to monitor
attach.jags(mod_lm_intercept)

mod_lm_intercept

# Now we can make plots of posterior values
par(mfrow = c(2,1))
hist(mu,40,col="grey",xlab="Mean",main="")
hist(sd.obs,40,col="grey",xlab=expression(sigma[obs]),main="")

# Set up the function
createMcmcList = function(jagsmodel) {
  McmcArray = as.array(jagsmodel$BUGSoutput$sims.array)
  McmcList = vector("list",length=dim(McmcArray)[2])
  for(i in 1:length(McmcList)) McmcList[[i]] = as.mcmc(McmcArray[,i,])
  McmcList = mcmc.list(McmcList)
  return(McmcList)
}
# Run the majority of the diagnostics that CODA() offers
gelmanDiags = gelman.diag(createMcmcList(mod_lm_intercept),multivariate=F)
autocorDiags = autocorr.diag(createMcmcList(mod_lm_intercept))
gewekeDiags = geweke.diag(createMcmcList(mod_lm_intercept))
heidelDiags = heidel.diag(createMcmcList(mod_lm_intercept))
 
######################################################################################
# 2. MODIFY THE ERRORS TO BE AUTOCORRELATED 
# no covariates, so intercept only. 
######################################################################################
jagsscript = cat("
model {	
   # priors on parameters
   mu ~ dnorm(0, 0.01); # This is normal with mean = 0, sd = 1/sqrt(0.01)
   tau.obs ~ dgamma(0.001,0.001); # This is inverse gamma
   sd.obs <- 1/sqrt(tau.obs); # sd is treated as derived parameter
   theta ~ dunif(-1,1);
   tau.cor <- tau.obs * (1-theta*theta); # Var = sigma2 * (1-rho^2)
   
   # Jags is not vectorized, so we have to loop over observations
   epsilon[1] <- Y[1] - mu;
   predY[1] <- mu;
   for(i in 2:N) {
      # for each observation, we'll assume the data to be normally distributed around
      # a common mean
      predY[i] <- mu + theta * epsilon[i-1]; 
      Y[i] ~ dnorm(predY[i], tau.cor);
      epsilon[i] <- (Y[i] - mu) - theta*epsilon[i-1];
      # The above 2 lines are inefficient, and could also be written as Y[i] ~ dnorm(mu, tau);
   }
}  

",file="lmcor_intercept.txt")

jags.data = list("Y"=Wind,"N"=N)
jags.params=c("sd.obs","predY","mu","theta")
model.loc=("lmcor_intercept.txt")
mod_lmcor_intercept = jags(jags.data, parameters.to.save=jags.params, model.file=model.loc, 
n.chains = 3, n.burnin=5000, n.thin=1, n.iter=10000, DIC=TRUE)   
 
 
######################################################################################
# 3. MAKE THE MODEL AN MA(1) MODEL WITH MA COEFFICIENT ESTIMATED
# no covariates. We can drop the global mean 'mu' and for clarity we'll rename tau.resid 
# to tau.obs, because the error is obs error variation. Note too that we define predY[1]
######################################################################################
jagsscript = cat("
model {	
   # priors on parameters
   mu ~ dnorm(0, 0.01); # This is normal with mean = 0, sd = 1/sqrt(0.01)
   tau.obs ~ dgamma(0.001,0.001); # This is inverse gamma
   sd.obs <- 1/sqrt(tau.obs); # sd is treated as derived parameter

   theta ~ dnorm(0,1);
   
   # Jags is not vectorized, so we have to loop over observations
   epsilon[1] <- Y[1] - mu; # error for first time point
   predY[1] <- mu;
   for(i in 2:N) {
      # for each observation, we'll assume the data to be normally distributed around
      # a common mean
      predY[i] <- mu + theta*epsilon[i-1];
      Y[i] ~ dnorm(predY[i], tau.obs);
      epsilon[i] <- (Y[i] - mu) - theta*epsilon[i-1];      
      # The above 2 lines are inefficient, and could also be written as Y[i] ~ dnorm(mu, tau);
   }
}  

",file="ma1_intercept.txt")

jags.data = list("Y"=Wind,"N"=N)
jags.params=c("sd.obs","predY","theta")
model.loc=("ma1_intercept.txt")
mod_ma1_intercept = jags(jags.data, parameters.to.save=jags.params, model.file=model.loc, n.chains = 3, n.burnin=5000, n.thin=1, n.iter=10000, DIC=TRUE)  
 
 
######################################################################################
# 4. MAKE THE MODEL AN AR(1) MODEL WITH NO ESTIMATED AR COEFFICIENT = RANDOM WALK
# no covariates. The model is y[t] ~ Normal(mu + y[n-1], sigma) for clarity we'll call the prcsn
# tau.pro, because the error is process variation. Note too that we have to define predY[1]
######################################################################################
jagsscript = cat("
model {	
   # priors on parameters
   mu ~ dnorm(0, 0.01); # This is normal with mean = 0, sd = 1/sqrt(0.01)
   tau.pro ~ dgamma(0.001,0.001); # This is inverse gamma
   sd.pro <- 1/sqrt(tau.pro); # sd is treated as derived parameter

   # Jags is not vectorized, so we have to loop over observations
   predY[1] <- mu;
   for(i in 2:N) {
      # for each observation, we'll assume the data to be normally distributed around
      # a common mean
      predY[i] <- mu + Y[i-1]; 
      Y[i] ~ dnorm(predY[i], tau.pro);
      # The above 2 lines are inefficient, and could also be written as Y[i] ~ dnorm(mu + Y[i-1], tau);
   }
}  

",file="rw_intercept.txt")

jags.data = list("Y"=Wind,"N"=N)
jags.params=c("sd.pro","predY","mu")
model.loc=("rw_intercept.txt")
mod_rw_intercept = jags(jags.data, parameters.to.save=jags.params, model.file=model.loc, n.chains = 3, n.burnin=5000, n.thin=1, n.iter=10000, DIC=TRUE)  
      
      
######################################################################################
# 5. MAKE THE MODEL AN AR(1) MODEL WITH AND ESTIMATED AR COEFFICIENT
# no covariates. We're introducting a new AR coefficient 'theta', rather than just modeling the process
# as a random walk, so the model is y[t] ~ Normal(mu + theta*y[n-1], sigma) for clarity we'll call the prcsn
# tau.pro, because the error is process variation. Note too that we have to define predY[1]
######################################################################################
jagsscript = cat("
model {	
   # priors on parameters
   mu ~ dnorm(0, 0.01); # This is normal with mean = 0, sd = 1/sqrt(0.01)
   tau.pro ~ dgamma(0.001,0.001); # This is inverse gamma
   sd.pro <- 1/sqrt(tau.pro); # sd is treated as derived parameter
   theta ~ dnorm(0, 1); # this is the ar coefficient
   
   # Jags is not vectorized, so we have to loop over observations
   predY[1] <- Y[1];
   for(i in 2:N) {
      # for each observation, we'll assume the data to be normally distributed around
      # a common mean
      predY[i] <- mu + theta * Y[i-1]; 
      Y[i] ~ dnorm(predY[i], tau.pro);
      # The above 2 lines are inefficient, and could also be written as Y[i] ~ dnorm(mu + theta * Y[i-1], tau);
   }
}  

",file="ar1_intercept.txt")

jags.data = list("Y"=Wind,"N"=N)
jags.params=c("sd.pro","predY","mu","theta")
model.loc=("ar1_intercept.txt")
mod_ar1_intercept = jags(jags.data, parameters.to.save=jags.params, model.file=model.loc, n.chains = 3, n.burnin=5000, n.thin=1, n.iter=10000, DIC=TRUE)  

######################################################################################
# 6. MAKE THE SS MODEL a univariate random walk
# no covariates. 
######################################################################################
jagsscript = cat("
model {	
   # priors on parameters
   mu ~ dnorm(0, 0.01); # This is normal with mean = 0, sd = 1/sqrt(0.01)
   tau.pro ~ dgamma(0.001,0.001); # This is inverse gamma
   sd.pro <- 1/sqrt(tau.pro); # sd is treated as derived parameter
   tau.obs ~ dgamma(0.001,0.001); # This is inverse gamma
   sd.obs <- 1/sqrt(tau.obs); # sd is treated as derived parameter
      
   # Jags is not vectorized, so we have to loop over observations
   X[1] ~ dnorm(0,0.001);
   predY[1] <- X[1];
   Y[1] ~ dnorm(X[1], tau.obs);

   for(i in 2:N) {
      predX[i] <- X[i-1] + mu; 
      X[i] ~ dnorm(predX[i],tau.pro); # Process variation
      predY[i] <- X[i];
      Y[i] ~ dnorm(X[i], tau.obs); # Observation variation
   }
}  

",file="ss_model.txt")

jags.data = list("Y"=Wind,"N"=N)
jags.params=c("sd.pro","sd.obs","predY","mu")
model.loc=("ss_model.txt")
mod_ss = jags(jags.data, parameters.to.save=jags.params, model.file=model.loc, n.chains = 3, n.burnin=5000, n.thin=1, n.iter=10000, DIC=TRUE)  
   
   
######################################################################################
# Optional: Let's go back to our original model and include some covariates
# We'll use temperature as a predictor of wind in the air quality dataset. 
######################################################################################
jagsscript = cat("
model {	
   # priors on parameters
   mu ~ dnorm(0, 0.01); # This is normal with mean = 0, sd = 1/sqrt(0.01)
   tau.obs ~ dgamma(0.001,0.001); # This is inverse gamma
   sd.obs <- 1/sqrt(tau.obs); # sd is treated as derived parameter
   B ~ dnorm(0,0.01);
   
   # Jags is not vectorized, so we have to loop over observations
   for(i in 1:N) {
      # for each observation, we'll assume the data to be normally distributed around
      # a common mean
      predY[i] <- mu + C[i]*B; 
      Y[i] ~ dnorm(predY[i], tau.obs);
      # The above 2 lines are inefficient, and could also be written as Y[i] ~ dnorm(mu, tau);
   }
}  

",file="lm.txt")

jags.data = list("Y"=Wind,"N"=N,"C"=Temp)
jags.params=c("sd.obs","predY","mu","B")
model.loc=("lm.txt")
mod_lm = jags(jags.data, parameters.to.save=jags.params, model.file=model.loc, 
n.chains = 3, n.burnin=5000, n.thin=1, n.iter=10000, DIC=TRUE)  
######################################################################################
# If you have extra time, go back and try to include Temperature in some of the other
# models
######################################################################################

# Linear regression with autocorrelated errors'
# Non-linear model, e.g. Ricker
# 
    
   
   
   
   
   
   
   
plotModelOutput = function(jagsmodel, Y) {
	# attach the model
	attach.jags(jagsmodel)
	x = seq(1,length(Y))
	summaryPredictions = cbind(apply(predY,2,quantile,0.025), apply(predY,2,mean), apply(predY,2,quantile,0.975))
	plot(Y, col="white",ylim=c(min(c(Y,summaryPredictions)),max(c(Y,summaryPredictions))),xlab="",ylab="95% CIs of predictions and data",main=paste("JAGS results:",jagsmodel$model.file))
	polygon(c(x,rev(x)), c(summaryPredictions[,1], rev(summaryPredictions[,3])), col="grey70",border=NA)
	lines(summaryPredictions[,2])
	points(Y)
}

createMcmcList = function(jagsmodel) {
	McmcArray = as.array(jagsmodel$BUGSoutput$sims.array)
	McmcList = vector("list",length=dim(McmcArray)[2])
	for(i in 1:length(McmcList)) McmcList[[i]] = as.mcmc(McmcArray[,i,])
	McmcList = mcmc.list(McmcList)
	return(McmcList)
}



print(plotModelOutput(mod_lm_intercept, Y=Wind))
print(plotModelOutput(mod_ma1_intercept, Y=Wind))
print(plotModelOutput(mod_ar1a_intercept, Y=Wind))
print(plotModelOutput(mod_ar1b_intercept, Y=Wind))


