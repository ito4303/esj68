functions {
  /**
   * Return log probability of N-mixture model for a site
   * 
   * @param count    Count
   * @param num_rep  Number of replications
   * @param max_n    Maximum abundance
   * @param lambda   Mean abundance
   * @param p        Detection probability
   *
   * @return         Log probability
   */
  real n_mixture_lpmf(int[] count, int max_n,
                      real lambda, real p) {
                 
    int c_max = max(count);
    vector[max_n + 1] lp;

    for (k in 0:(c_max - 1))
      lp[k + 1] = negative_infinity();
    for (k in c_max:max_n) 
      lp[k + 1] = poisson_lpmf(k | lambda) + binomial_lpmf(count | k, p);
    return log_sum_exp(lp);
  }

  real partial_sum(int[] site,
                   int start, int end,
                   int[, ] count,
                   int max_n, real lambda, real p) {
    real lp = 0;

    for (m in start:end)
      lp = lp + n_mixture_lpmf(count[m] | max_n, lambda, p);
    return lp;
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
  int site[M] = rep_array(0, M); //dummy for site index
}

parameters {
  real<lower = 0> lambda;       // mean abundance
  real<lower = 0, upper = 1> p; // detection probability
}

model {
  int grainsize = 1;

  lambda ~ normal(0, 20);
  target += reduce_sum(partial_sum, site, grainsize, C, K, lambda, p);
}
