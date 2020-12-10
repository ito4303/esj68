data {
  int<lower = 0> M;
  int<lower = 0> J;
  int<lower = 0, upper = 1> Y[M, J];
}

transformed data {
  int<lower = 0> Ysum[M];

  for (m in 1:M)
    Ysum[m] = sum(Y[m, ]);
}

parameters {
  real<lower = 0, upper = 1> psi;
  real<lower = 0, upper = 1> p;
}

model {
  for (m in 1:M) {
    if (Ysum[m] > 0) { // detected
      target += bernoulli_lpmf(1 | psi)
              + bernoulli_lpmf(Y[m, ] | p);
    } else {  // not detected
      real lp[2];

      lp[1] = bernoulli_lpmf(0 | psi);
      lp[2] = bernoulli_lpmf(1 | psi)
            + bernoulli_lpmf(Y[m, ] | p);
      target += log_sum_exp(lp);
    }
  }
}
