data {
  int<lower = 0> M;       // number of sites
  int<lower = 0> J;       // number of replications
  int<lower = 0> K;       // upper limit of the abundance
  int<lower = 0> N[M];    // abundance for each site
  int<lower = 0> C[M, J]; // count for each site and replication
}

parameters {
  real<lower = 0> lambda;       // mean abundance
  real<lower = 0, upper = 1> p; // detection probability
}

model {
  for (m in 1:M) {
    int c_max = max(C[m]);
    vector[K + 1] lp;

    for (k in 0:(c_max - 1))
      lp[k + 1] = negative_infinity();
    for (k in c_max:K) 
      lp[k + 1] = poisson_lpmf(k | lambda) + binomial_lpmf(C[m] | k, p);
    target += log_sum_exp(lp);
  }
}
