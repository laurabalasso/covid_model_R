library(rstan)
library(bayesplot)
source('R/data_cleaning_it.R')
source('R/regional_model_data.R')
source('R/rt_plots.R')

### Model for Lombardia


data_lombardia <- get_model_data(data_it, 'Lombardia', initial_date = as.Date('2020-08-30'))

p_delay <- get_delay_distribution()

nonzero_days_l <- which(data_lombardia$total != 0)

stan_data_l <- list(N = nrow(data_lombardia),
                  conv_gt = get_gt_convolution_ln2(nrow(data_lombardia)),
                  length_delay = length(p_delay),
                  p_delay = p_delay,
                  exposures = exposures_from_total(data_lombardia$total),
                  N_nonzero = length(nonzero_days_l),
                  nonzero_positives = data_lombardia$positive[nonzero_days_l],
                  nonzero_days = nonzero_days_l
)


compiled_model <- stan_model('stan/rt_model.stan')

fit_model_lomb <- sampling(compiled_model, data=stan_data_l, iter = 200)


print(fit_model_lomb, pars = 'r_t')

## trace plots

mcmc_trace(as.array(fit_model_lomb, pars = c('r_t[1]', 'r_t[10]', 'r_t[50]', 'r_t[60]')),
           np = nuts_params(fit_model_lomb)
           )

mcmc_trace(as.array(fit_model_lomb, pars = c('eta[1]', 'eta[10]', 'eta[50]', 'eta[60]')),
           np = nuts_params(fit_model_lomb)
)




## effective sample size ratios 

ratios1 <- neff_ratio(fit_model_lomb, pars = c('eta'))
mcmc_neff(ratios1) + yaxis_text(hjust = 1)

ratios2 <- neff_ratio(fit_model_lomb, pars = c('r_t'))
mcmc_neff(ratios2) + yaxis_text(hjust = 1)

check_n_eff(fit_model_lomb)



## posterior predictive check

y_rep <- as.matrix(fit_model_lomb, pars = "y_rep")
ppc_dens_overlay(y = data_lombardia$positive[nonzero_days_l], y_rep[1:200, ]) 

ppc_intervals(
  y = data_lombardia$positive[nonzero_days_l],
  yrep = y_rep
)




### R_t curve


plot_rt(data_lombardia, fit_model_lomb)


