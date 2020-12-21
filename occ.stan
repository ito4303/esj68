data {
  int<lower = 0> M;
  int<lower = 0> J;
  int<lower = 0, upper = 1> Y[M, J];
  vector[M] Elev;
  vector[M] Forest;
  matrix[M, J] Wind;
}

transformed data {
  int<lower = 0> Ysum[M];

  for (m in 1:M)
    Ysum[m] = sum(Y[m, ]);
}

parameters {
  vector[3] beta;
  vector[2] alpha;
}

transformed parameters {
  vector[M] logit_psi = beta[1]
                        + beta[2] * Elev
                        + beta[3] * Forest;
  matrix[M, J] logit_p = alpha[1]
                         + alpha[2] * Wind;
}

model {
  for (m in 1:M) {
    if (Ysum[m] > 0) { // detected
      target += bernoulli_logit_lpmf(1 | logit_psi[m])
              + bernoulli_logit_lpmf(Y[m, ] | logit_p[m, ]);
    } else {          // not detected
      real lp[2];

      lp[1] = bernoulli_logit_lpmf(0 | logit_psi[m]);
      lp[2] = bernoulli_logit_lpmf(1 | logit_psi[m])
            + bernoulli_logit_lpmf(Y[m, ] | logit_p[m, ]);
      target += log_sum_exp(lp);
    }
  }
}
