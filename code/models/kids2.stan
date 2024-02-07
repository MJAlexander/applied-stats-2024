data {
  int<lower=0> N;          // number of kids
  vector[N] y;             // scores
  real mu0;                // mean of mu prior
  real<lower=0> sigma0;   //  sd of mu prior
}
parameters {
  real mu; 
  real<lower=0> sigma;
}
transformed parameters {
}
model {
  //priors
  mu ~ normal(mu0, sigma0);
  sigma ~ normal(0,10);
  
  //target += normal_lpdf(y | mu, sigma);
  //equivalent:
  y ~ normal(mu, sigma);
}