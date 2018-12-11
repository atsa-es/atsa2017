data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real x0;
  real<lower=-0.999,upper=0.999> phi;
  real mu;
  vector[N-1] pro_dev;
  real<lower=0> sigma_process;
  real<lower=0> sigma_obs;
}
transformed parameters {
  vector[N] pred;
  pred[1] = x0;
  for(i in 2:N) {
    pred[i] = mu + phi*(pred[i-1] - mu) + pro_dev[i-1];
  }
}
model {
  x0 ~ normal(0,10);
  phi ~ normal(0,10);
  mu ~ normal(0,10);
  sigma_process ~ cauchy(0,5);
  sigma_obs ~ cauchy(0,5);
  pro_dev ~ normal(0, sigma_process);
  y ~ normal(pred, sigma_obs);
}