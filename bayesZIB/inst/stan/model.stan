functions {
  real zib_lpmf(int y, real w, real p) {
    return y*log(w*p)+(1-y)*log(1-w*p);
  }
}

data {
  // Define variables in data
  // Number of observations (an integer)
  int<lower=0> N;
  // Number of parameters (zero-inflated part)
  int<lower=0> M1;
  // Number of parameters (non zero-inflated part)
  int<lower=0> M2;
  // Variables
  int<lower=0, upper=1> y[N]; // outcome
  matrix[N, M1] X1;           // design matrix (zero-inflated part)
  matrix[N, M2] X2;           // design matrix (non zero-inflated part)
  // Hyperparameters
  real<lower=0> s_theta[M1]; 
  real<lower=0> s[M2]; 
}

parameters {
  // define parameters to estimate
  vector[M1] theta;
  vector[M2] beta;
}

transformed parameters {
  // Probability trasformation from linear predictor
  vector[N] eta_theta;
  vector[N] eta_beta;
  real<lower=0.0, upper=1.0> w[N];
  real<lower=0.0, upper=1.0> p[N];
  
  eta_theta = X1 * theta;
  eta_beta  = X2 * beta;
 
  for (n in 1:N) {
    w[n] = exp(eta_theta[n])/(1+exp(eta_theta[n]));
    p[n] = exp(eta_beta[n])/(1+exp(eta_beta[n]));
  }
}

model {
  // Priors
  for (i in 1:M1) {
    target += normal_lpdf(theta[i] | 0, s_theta[i]);
  }
  for (i in 1:M2) {
   target += normal_lpdf(beta[i] | 0, s[i]);
  }
  // Likelihood
  for (n in 1:N) {
    target += zib_lpmf(y[n] | w[n], p[n]);
  }
}
