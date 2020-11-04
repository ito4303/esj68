functions {
  /**
   * Returns log likelihood of N-mixture model
   * with 2 replicated observations using
   * bivariate Poisson distibution
   *
   * References
   * Dennis et al. (2015) Computational aspects of N-mixture models.
   *   Biometrics 71:237--246. DOI:10.1111/biom.12246
   * Stan users mailing list
   *   https://groups.google.com/forum/#!topic/stan-users/9mMsp1oB69g
   *
   * @param n          Number of observed individuals
   * @param log_lambda Log of Poisson mean of population size
   * @param p          Detection probability
   *
   * return Log probability
  */
  real bivariate_poisson_log_lpmf(int[] n, real log_lambda, real p) {
    real s[min(n) + 1];
    real log_theta_1 = log_lambda + log(p) + log1m(p);
    real log_theta_0 = log_lambda + log(p) * 2;

    if (size(n) != 2)
      reject("Size of n must be 2.");
    if (p < 0 || p > 1)
      reject("p must be in [0,1].");
    for (u in 0:min(n))
      s[u + 1] = poisson_log_lpmf(n[1] - u | log_theta_1)
               + poisson_log_lpmf(n[2] - u | log_theta_1)
               + poisson_log_lpmf(u | log_theta_0);
    return log_sum_exp(s);
  }

  real partial_sum(int[] site,
                   int start, int end,
                   int[, ] count,
                   int max_n, real log_lambda, real p) {
    real lp = 0;

    for (m in start:end)
      lp = lp + bivariate_poisson_log_lpmf(count[m] | log_lambda, p);
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

transformed parameters {
  real log_lambda = log(lambda);
}

model {
  int grainsize = 1;

  lambda ~ normal(0, 20);
  target += reduce_sum(partial_sum, site, grainsize,
                       C, K, log_lambda, p);
}
