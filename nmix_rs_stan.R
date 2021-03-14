library(cmdstanr)
options(mc.cores = parallel::detectCores())

stan_data <- list(M = nmix_data$nsites, J = nmix_data$nvisits,
                  C = nmix_data$C,
                  Cov_abn = nmix_data$site.cov[, 2],
                  Cov_det = nmix_data$survey.cov,
                  K = max(nmix_data$C) + 100)
model_file <- "nmix_rs.stan"
mod_nmix_stan <- cmdstan_model(model_file)
fit_nmix_stan <- mod_nmix_stan$sample(
  data = stan_data, seed = 1,
  chains = 4, parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000,
  refresh = 500)
fit_nmix_stan$print(c("beta", "alpha", "Ntotal"))
