data {
  int<lower = 0> M;
  int<lower = 0> J;
  int<lower = 0> C[M, J];
  vector[M] Cov_abn;
  matrix[M, J] Cov_det;
  int<lower = 0> K;
}

parameters {
  vector[2] beta;
  vector[2] alpha;
}

transformed parameters {
  vector[M] log_lambda = beta[1] + beta[2] * Cov_abn;
  matrix[M, J] logit_p = alpha[1] + alpha[2] * Cov_det;
  matrix[M, K + 1] lp;

  for (m in 1:M) {
    int c_max = max(C[m]);

    for (k in 0:(c_max - 1))
      lp[m, k + 1] = negative_infinity();
    for (k in c_max:K)
      lp[m, k + 1] = poisson_log_lpmf(k | log_lambda[m])
                    + binomial_logit_lpmf(C[m, ] | k, logit_p[m]);
  }
}

model {
  for (m in 1:M) {
    target += log_sum_exp(lp[m]);
  }
  beta ~ normal(0, 10);
  alpha ~ normal(0, 10);
}

generated quantities {
  int<lower = 0> N[M];
  int<lower = 0> Ntotal;

  for (m in 1:M)
    N[m] = categorical_rng(softmax(lp[m]')) - 1;
  Ntotal = sum(N);
}
