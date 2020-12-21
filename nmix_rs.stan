functions {
  /**
   * Return log probability of N-mixture model for a site
   * 
   * @param count      Count
   * @param max_n      Maximum abundance
   * @param log_lambda Log of mean abundance
   * @param logit_p    Logit of detection probability
   *
   * @return         Log probability
   */
  real n_mixture_log_lpmf(int[] count, int max_n,
                          real log_lambda, vector logit_p) {
                 
    int c_max = max(count);
    vector[max_n + 1] lp;

    for (k in 0:(c_max - 1))
      lp[k + 1] = negative_infinity();
    for (k in c_max:max_n) 
      lp[k + 1] = poisson_log_lpmf(k | log_lambda)
                  + binomial_logit_lpmf(count | k, logit_p);
    return log_sum_exp(lp);
  }

  real partial_sum(int[] site,
                   int start, int end,
                   int[, ] count,
                   int max_n,
                   vector log_lambda, matrix logit_p) {
    real lp = 0;

    for (m in start:end)
      lp = lp + n_mixture_log_lpmf(count[m] |
                                   max_n, log_lambda[m], logit_p[m, ]');
    return lp;
  }
}

data {
  int<lower = 0> M;
  int<lower = 0> J;
  int<lower = 0> Y[M, J];
  vector[M] Cov_abn;
  matrix[M, J] Cov_det;
  int<lower = 0> Max_N;
}

transformed data {
  int site[M] = rep_array(0, M); // dummy site index
}

parameters {
  vector[2] beta;
  vector[2] alpha;
}

transformed parameters {
  vector[M] log_lambda = beta[1]
                         + beta[2] * Cov_abn;
  matrix[M, J] logit_p = alpha[1]
                         + alpha[2] * Cov_det;
}

model {
  int grainsize = 1;

  target += reduce_sum(partial_sum, site, grainsize, Y, Max_N,
                       log_lambda, logit_p);
  beta ~ normal(0, 10);
  alpha ~ normal(0, 10);
}
