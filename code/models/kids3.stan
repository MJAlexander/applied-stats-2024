data {
  int<lower=0> N;          // number of kids
  int<lower=0> K;          // number of covariates
  vector[N] y;             // scores
  matrix[N, K] X;           // design matrix
}
parameters {
  real alpha;
  vector[K] beta; 
  real<lower=0> sigma;
}
transformed parameters {
}
model {
  //priors
  alpha ~ normal(0, 100);
  beta ~ normal(0, 10);
  sigma ~ normal(0,10);
  
  //likelihood
  y ~ normal(alpha + X*beta, sigma);
}