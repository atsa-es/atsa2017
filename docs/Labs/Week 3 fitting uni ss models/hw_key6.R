library(forecast)
#Question 6. Finally, we'll evaluate the predictive accuracy of forecasts using the forecast package using the 'airmiles' dataset,
#Load the data to use as follows:

dat=log(airmiles)
n = length(dat)
training.dat = dat[1:(length(airmiles)-3)]
test.dat = dat[(length(airmiles)-2):n]

#We'll use the last 3 data points to validate our fits. We want you to fit the following four models, and calculate the MASE statistic for each (using the forecasts for the 3 data points): ARIMA(0,0,0), ARIMA(1,0,0), ARIMA(0,0,1), ARIMA(1,0,1). Present them in a table similar to what Eric showed in lecture. Hint: use the forecast function and accuracy function to calculate the MASE statistic. 
accuracy(forecast(Arima(training.dat, order =c(0,0,0)), h=3), test.dat) # MASE = 9.129155
accuracy(forecast(Arima(training.dat, order =c(1,0,0)), h=3), test.dat) # MASE = 0.6906049
accuracy(forecast(Arima(training.dat, order =c(0,0,1)), h=3), test.dat) # MASE = 7.735312
accuracy(forecast(Arima(training.dat, order =c(1,0,1)), h=3), test.dat) # MASE = 0.5119703

# What this shows is that the ARMA(1,1) is the best, and the AR component strongly improves predictions
