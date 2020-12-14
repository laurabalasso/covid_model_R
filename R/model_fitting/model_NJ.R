library(rstan)
library(bayesplot)
source('R/data_cleaning_us.R')
source('R/regional_model_data.R')
source('R/rt_plots.R')

### model for New Jersey 

NJ_data <- get_model_data(data, 'NJ', as.Date("2020-06-30") )

p_delay <- get_delay_distribution()

nonzero_days_NJ <- which(NJ_data$total != 0)

stan_data_NJ <- list(N = nrow(NJ_data),
                  conv_gt = get_gt_convolution_ln2(nrow(NJ_data)),
                  length_delay = length(p_delay),
                  p_delay = p_delay,
                  exposures = exposures_from_total(NJ_data$total),
                  N_nonzero = length(nonzero_days_NJ),
                  nonzero_positives = NJ_data$positive[nonzero_days_NJ],
                  nonzero_days = nonzero_days_NJ
                  
)


compiled_model <- stan_model('stan/rt_model.stan')

fit_model_NJ <- sampling(compiled_model, data=stan_data_NJ, iter = 2000)

y_rep <- as.matrix(fit_model_NJ, pars = "y_rep")
ppc_dens_overlay(y = stan_data_NJ$nonzero_positives,
                 yrep = y_rep[1:1000, ])

ppc_intervals(y = stan_data_NJ$nonzero_positives,
              yrep = y_rep[1:1000, ])

fit_summary_NJ <- summary(fit_model_NJ)

## rt plot

plor_rt(NJ_data, fit_model_NJ,'NJ')
