library(rstan)
library(bayesplot)
library(ggplot2)
source('data_cleaning_it.R')
source('utility_functions.R')

### Model for Lombardia


data_lombardia <- get_model_data(data_it, 'Lombardia', initial_date = as.Date('2020-06-30'))

p_delay <- get_delay_distribution()

nonzero_days_l <- which(data_lombardia$total != 0)

stan_data_l <- list(N = nrow(data_lombardia),
                  conv_gt = get_gt_convolution(nrow(data_lombardia)),
                  length_delay = length(p_delay),
                  p_delay = p_delay,
                  exposures = exposures_from_total(data_lombardia$total),
                  N_nonzero = length(nonzero_days_l),
                  nonzero_positives = data_lombardia$positive[nonzero_days_l],
                  nonzero_days = nonzero_days_l
)


compiled_model <- stan_model('rt_model.stan')

fit_model_lomb <- sampling(compiled_model, data=stan_data_l, iter = 2000)


print(fit_model_lomb, pars = 'r_t')
print(fit_model_lomb, pars = 'mu')

## trace plots

mcmc_trace(as.array(fit_model_lomb, pars = c('r_t[1]', 'r_t[10]', 'r_t[50]', 'r_t[100]')),
           np = nuts_params(fit_model_lomb)
           )

mcmc_trace(as.array(fit_model_lomb, pars = c('mu[1]', 'mu[10]', 'mu[50]', 'mu[100]')),
           np = nuts_params(fit_model_lomb)
)


mcmc_trace(
  as.array(fit_model_lomb, pars = c('alpha', 'seed')),
  np = nuts_params(fit_model_lomb)
)



## effective sample size ratios 

ratios1 <- neff_ratio(fit_model_lomb, pars = c('mu'))
mcmc_neff(ratios1) + yaxis_text(hjust = 1)

ratios2 <- neff_ratio(fit_model_lomb, pars = c('r_t'))
mcmc_neff(ratios2) + yaxis_text(hjust = 1)

check_n_eff(fit_model_lomb)



## posterior predictive check

y_rep <- as.matrix(fit_model_lomb, pars = "y_rep")
ppc_dens_overlay(y = data_lombardia$positive[nonzero_days_l], y_rep[1:1000, ]) 

ppc_intervals(
  y = data_lombardia$positive[nonzero_days_l],
  yrep = y_rep,
  x = data_lombardia$date[nonzero_days_l]
)


### R_t curve

fit_summary_lomb <- summary(fit_model_lomb)

medians_rt <- fit_summary_lomb$summary[, '50%'][140: (140+136)]
min_rt_50_interval <- fit_summary_lomb$summary[, '25%'][140: (140+136)]
max_rt_50_interval <- fit_summary_lomb$summary[, '75%'][140: (140+136)]
min_rt_95_interval <- fit_summary_lomb$summary[, '2.5%'][140: (140+136)]
max_rt_95_interval <- fit_summary_lomb$summary[, '97.5%'][140: (140+136)]


ggplot(data = NULL, aes(x = data_lombardia$date, y = medians_rt)) + 
  geom_line() + 
  xlab('Date') +
  ylab('') +
  ggtitle( 'Lombardia r_t')+
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  geom_vline(xintercept = data_lombardia$date[1]) +
  geom_ribbon(aes(ymin = min_rt_50_interval, ymax = max_rt_50_interval), alpha= 0.5, fill = 'darkred') +
  geom_ribbon(aes(ymin = min_rt_95_interval, ymax = max_rt_95_interval), alpha= 0.1, fill = 'darkred')




