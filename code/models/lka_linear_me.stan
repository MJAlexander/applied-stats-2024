
data {
  int<lower=0> N; // number of observations
  int<lower=0> T; //number of years
  int<lower=0> mid_year; // mid-year of study
  vector[N] y; //log ratio
  vector[N] se; // standard error around observations
  vector[T] years; // unique years of study
  int<lower=0> year_i[N]; // year index of observations
  
}

parameters {
  real alpha;
  real beta;

}

transformed parameters{
  vector[T] mu;
  
  for(t in 1:T){
    mu[t] = alpha + beta*(years[t] - mid_year);
  }
}

model {
  
  y ~ normal(mu[year_i], se);
  
  alpha ~ normal(0, 1);
  beta ~ normal(0,1);
}

