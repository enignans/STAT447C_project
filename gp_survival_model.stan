data {
  int<lower=1> N;                
  int<lower=1> D;                 
  matrix[N, D] X;                 
  vector[N] y;                    
  int<lower=0,upper=1> status[N]; 
}

parameters {
  real<lower=0> alpha;            
  vector<lower=0>[D] ell;         
  real<lower=0> sigma;            
  vector[N] f_raw;               
}

transformed parameters {
  matrix[N, N] K;
  matrix[N, N] L_K;
  vector[N] f;

  for (i in 1:N) {
    for (j in i:N) {
      real sqdist = 0;
      for (d in 1:D) {
        sqdist += square((X[i,d] - X[j,d]) / ell[d]);
      }
      K[i,j] = square(alpha) * exp(-0.5 * sqdist);
      if (i != j) {
        K[j,i] = K[i,j];
      }
    }
    K[i,i] += 1e-6;  
  }

  L_K = cholesky_decompose(K);
  f = L_K * f_raw;
}

model {
  // Priors
  alpha ~ normal(0, 1);       
  ell ~ lognormal(0, 1);    
  sigma ~ normal(0, 1);       

  // GP latent function prior
  f_raw ~ normal(0, 1);

  // Likelihood
  for (n in 1:N) {
    if (status[n] == 1) {
      y[n] ~ normal(f[n], sigma); 
    } else {
      target += normal_lccdf(y[n] | f[n], sigma); 
    }
  }
}

generated quantities {
  vector[N] y_rep;
  for (n in 1:N) {
    y_rep[n] = normal_rng(f[n], sigma);
  }
}
