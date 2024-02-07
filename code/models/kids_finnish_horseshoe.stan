data {
  int<lower=0> N;          // number of kids
  int<lower=0> K;          // number of covariates
  vector[N] y;             // scores
  matrix[N, K] X;           // design matrix
  real <lower=0> slab_scale ; // slab scale for the regularized horseshoe
  real <lower=0> slab_df ; // slab degrees of freedom for the regularized horseshoe
}
parameters {
  real alpha;
  vector[K] z; 
  real<lower=0> sigma;
  real<lower=0> tau;
  vector<lower=0>[K] lambda;
  real<lower=0> caux;
}
transformed parameters {
  vector<lower=0>[K] lambda_tilde;
  vector[K] beta;
  vector[N] f;
  real<lower=0> c;
  c = slab_scale * sqrt ( caux );
  lambda_tilde = sqrt ( c ^2 * square ( lambda ) ./ (c ^2 + tau ^2* square ( lambda )) );
  beta = z .* lambda_tilde*tau;
  f = alpha + X*beta;
}

model {
  //priors
  alpha ~ normal(0, 100);
  z ~ normal(0,1);
  sigma ~ normal(0,10);
  lambda ~ student_t(2, 0, 1);
  tau ~ student_t(2, 0, 2*sigma);
  caux ~ inv_gamma(0.5* slab_df , 0.5* slab_df );
  
  //likelihood
  y ~ normal(f, sigma);
}