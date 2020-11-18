library(unmarked)
library(cmdstanr)
options(mc.cores = parallel::detectCores())
set_cmdstan_path("~/.cmdstan/cmdstan-2.25.0")

lambda <- 2.5 # mean abundance
p <- 0.4      # detection probability
M <- 150      # Number of sites
J <- 2        # Number of replications

set.seed(123)
N <- rpois(M, lambda) # abundance for each site
C <- matrix(0, nrow = M, ncol = J) # count for each site
for (j in 1:J)
  C[, j] <-rbinom(M, N, p)

## Stan

### Simple implimentation

stan_data <- list(M = M, J = J, K = 100, N = N, C = C)
stan_model1 <- cmdstan_model("simplest_N-mixture1.stan")
fit_stan1 <- stan_model1$sample(data = stan_data,
                                chains = 4, parallel_chains = getOption("mc.cores", 1),
                                iter_warmup = 1000, iter_sampling = 1000,
                                refresh = 500)
fit_stan1$summary(variables = c("lambda", "p"))
fit_stan1$time()

### use reduce_sum

stan_data <- list(M = M, J = J, K = 100, N = N, C = C)
stan_model2 <- cmdstan_model("simplest_N-mixture2.stan",
                             cpp_options = list(stan_threads = TRUE))
fit_stan2 <- stan_model2$sample(data = stan_data,
                                chains = 4,
                                parallel_chains = getOption("mc.cores", 1),
                                threads_per_chain = 3,
                                iter_warmup = 1000, iter_sampling = 1000,
                                refresh = 500)
fit_stan2$cmdstan_diagnose()
fit_stan2$summary(variables = c("lambda", "p"))
fit_stan2$time()
