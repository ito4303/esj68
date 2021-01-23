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
  real bivariate_poisson_log_lpmf(int[] n, real log_lambda, vector p) {
    real s[min(n) + 1];
    real log_theta_2 = log_lambda + log(p[1]) + log1m(p[2]);
    real log_theta_1 = log_lambda + log1m(p[1]) + log(p[2]);
    real log_theta_0 = log_lambda + log(p[1]) + log(p[2]);

    if (size(n) != 2 || num_elements(p) != 2)
      reject("Size of n and p must be 2.");
    if (p[1] < 0 || p[1] > 1 || p[2] < 0 || p[2] > 1)
      reject("p must be in [0,1].");
    for (u in 0:min(n))
      s[u + 1] = poisson_log_lpmf(n[1] - u | log_theta_2)
               + poisson_log_lpmf(n[2] - u | log_theta_1)
               + poisson_log_lpmf(u | log_theta_0);
    return log_sum_exp(s);
  }

  int n_mixture_rng(int[] count, int max_n,
                    real log_lambda, vector p) {
    int c_max = max(count);
    vector[max_n + 1] lp;

    for (k in 0:(c_max - 1))
      lp[k + 1] = negative_infinity();
    for (k in c_max:max_n) 
      lp[k + 1] = poisson_log_lpmf(k | log_lambda)
                  + binomial_lpmf(count | k, p);
    return categorical_rng(softmax(lp)) - 1;
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

parameters {
  vector[2] beta;
  vector[2] alpha;
}

transformed parameters {
  vector[M] log_lambda = beta[1] + beta[2] * Cov_abn;
  matrix[M, J] p = inv_logit(alpha[1] + alpha[2] * Cov_det);
}

model {
  for (m in 1:M)
    Y[m, ] ~ bivariate_poisson_log(log_lambda[m], p[m]');
  beta ~ normal(0, 10);
  alpha ~ normal(0, 10);
}

generated quantities {
  int<lower = 0> N[M];
  int<lower = 0> Ntotal;

  for (m in 1:M)
    N[m] = n_mixture_rng(Y[m, ], Max_N, log_lambda[m], p[m]');
  Ntotal = sum(N);
}
