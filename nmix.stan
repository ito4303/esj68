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
  vector[M] log_lambda = beta[1]
                         + beta[2] * Cov_abn;
  matrix[M, J] logit_p = alpha[1]
                         + alpha[2] * Cov_det;
}

model {
  for (m in 1:M) {
    int y_max = max(Y[m]);
    vector[Max_N + 1] lp;

    for (n in 0:(y_max - 1))
      lp[n + 1] = negative_infinity();
    for (n in y_max:Max_N)
      lp[n + 1] = poisson_log_lpmf(n | log_lambda[m])
                  + binomial_logit_lpmf(Y[m, ] | n, logit_p[m, ]);
    target += log_sum_exp(lp);
  }
  beta ~ normal(0, 10);
  alpha ~ normal(0, 10);
}
