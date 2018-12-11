require(coda)
require(R2jags)

#set-up the data
data(airquality)
Wind = airquality$Wind # wind speed
N = dim(airquality)[1] # number of data points

#Define the data and parameters to track
jags.data = list("Y"=Wind,"N"=N) # named list of data
jags.params=c("sd.obs","mu") # parameters in the linear regression model

# LINEAR REGRESSION Intercept only
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

## Fit the model
mod_lm_intercept = jags(jags.data, parameters.to.save=jags.params, 
                        model.file=model.loc, n.chains = 3, n.burnin=5000, n.thin=1, n.iter=10000, DIC=TRUE)  

# Now we can make plots of posterior values
attach.jags(mod_lm_intercept)
par(mfrow = c(2,1))
hist(mu,40,col="grey",xlab="Mean",main="")
hist(sd.obs,40,col="grey",xlab=expression(sigma[obs]),main="")



