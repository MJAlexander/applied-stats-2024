
data {
  int<lower=0> N; //number of years
  int<lower=0> S; //number of states
  int<lower=0> K; //knot points
  matrix[N,S] y; //log of entries per capita
  matrix[N,K] B; //B splines
}


parameters {
  matrix[K,S] alpha;
  vector<lower=0>[S] sigma;
  vector<lower=0>[S] sigma_alpha;
  real mu_sigma;
  real<lower=0> tau;
}

transformed parameters {
  matrix[N,S] mu;
  vector[S] log_sigma_alpha;
  
  for(i in 1:N){
    for(s in 1:S){
      mu[i,s] = B[i,]*alpha[,s];
    }
  }
  
  for(s in 1:S){
    log_sigma_alpha[s] = log(sigma_alpha[s]);
  }
  
}

model {
  for(s in 1:S){
   y[,s] ~ normal(mu[,s], sigma[s]); 
  alpha[1,s] ~ normal(0, sigma_alpha[s]);
  alpha[2,s] ~ normal(alpha[1,s], sigma_alpha[s]);
  alpha[3:K,s] ~ normal(2*alpha[2:(K-1),s] - alpha[1:(K-2),s], sigma_alpha[s]);
  sigma[s] ~ normal(0,1);
  log_sigma_alpha[s] ~ normal(mu_sigma, tau);
  mu_sigma ~ normal(0,1);
  tau ~ normal(0,1);
  }
}

