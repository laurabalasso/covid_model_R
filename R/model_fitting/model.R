library(rstan)
library(bayesplot)
source('R/data_cleaning_us.R')
source('R/regional_model_data.R')
source('R/rt_plots.R')


### model for Oregon from March to June 

X <- get_model_data(data, 'OR', as.Date('2020-06-30'))

p_delay <- get_delay_distribution()

nonzero_days <- which(X$total != 0)

stan_data <- list(N = nrow(X),
                  conv_gt = get_gt_convolution_ln2(nrow(X)),
                  length_delay = length(p_delay),
                  p_delay = p_delay,
                  exposures = exposures_from_total(X$total),
                  N_nonzero = length(nonzero_days),
                  nonzero_positives = X$positive[nonzero_days],
                  nonzero_days = nonzero_days
                  
)


compiled_model <- stan_model('stan/rt_model.stan')

fit_model <- sampling(compiled_model, data=stan_data, iter = 2000)


print(fit_model, pars = 'r_t')
print(fit_model, pars = 'mu')


## trace plots

mcmc_trace(
  as.array(fit_model,pars = c('r_t[1]', 'r_t[10]', 'r_t[50]', 'r_t[100]')),
  np = nuts_params(fit_model)
)

mcmc_trace(
  as.array(fit_model,pars = c('mu[1]','mu[10]', 'mu[50]', 'mu[100]' )),
  np = nuts_params(fit_model)
)

mcmc_trace(
  as.array(fit_model, pars = c('alpha', 'seed')),
  np = nuts_params(fit_model)
)


## effective sample size ratios 

ratios1 <- neff_ratio(fit_model, pars = c('mu'))
mcmc_neff(ratios1) + yaxis_text(hjust = 1)

ratios2 <- neff_ratio(fit_model, pars = c('r_t'))
mcmc_neff(ratios2) + yaxis_text(hjust = 1)

## posterior predictive check

y_rep <- as.matrix(fit_model, pars = "y_rep")
ppc_dens_overlay(y = X$positive[nonzero_days], y_rep[1:1000, ])


## r_t curve

fit_summary <- summary(fit_model)

plot_rt(X, fit_model, 'Oregon')
