functions {
  real partial_sum(int[] dummy_k,
                   int start, int end,
                   int[] count,
                   real lambda,
                   real p) {
    vector[end - start + 1] lp;

    for (k in start:end) {
      if (k - 1 < max(count))
        lp[k - start + 1] = negative_infinity();
      else
        lp[k - start + 1] = poisson_lpmf(k - 1 | lambda)
                            + binomial_lpmf(count | k - 1, p);
    }
    return log_sum_exp(lp);
  }
}

data {
  int<lower = 0> M;       // number of sites
  int<lower = 0> J;       // number of replications
  int<lower = 0> K;       // upper limit of the abundance
  int<lower = 0> N[M];    // abundance for each site
  int<lower = 0> C[M, J]; // count for each site and replication
}

transformed data {
  int dummy_k[K + 1] = rep_array(0, K + 1);
}

parameters {
  real<lower = 0> lambda;       // mean abundance
  real<lower = 0, upper = 1> p; // detection probability
}

model {
  int grainsize = 50;

  lambda ~ normal(0, 5);
  for (m in 1:M) {
    target += reduce_sum_static(partial_sum, dummy_k, grainsize, C[m], lambda, p);
//    target += partial_sum(dummy_k, 1, K + 1, C[m], lambda, p);

  }
}
