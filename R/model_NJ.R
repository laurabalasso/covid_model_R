library(rstan)
library(bayesplot)
library(ggplot2)
source('data_cleaning_us.R')
source('utility_functions.R')

### model for New Jersey 

NJ_data <- get_model_data(data, 'NJ', as.Date("2020-06-30") )

p_delay <- get_delay_distribution()

nonzero_days_NJ <- which(NJ_data$total != 0)

stan_data_NJ <- list(N = nrow(NJ_data),
                  conv_gt = get_gt_convolution(nrow(NJ_data)),
                  length_delay = length(p_delay),
                  p_delay = p_delay,
                  exposures = exposures_from_total(NJ_data$total),
                  N_nonzero = length(nonzero_days_NJ),
                  nonzero_positives = NJ_data$positive[nonzero_days_NJ],
                  nonzero_days = nonzero_days_NJ
                  
)


compiled_model <- stan_model('rt_model.stan')

fit_model_NJ <- sampling(compiled_model, data=stan_data_NJ, iter = 2000)

y_rep <- as.matrix(fit_model_NJ, pars = "y_rep")
ppc_dens_overlay(y = stan_data_NJ$nonzero_positives,
                 yrep = y_rep[1:1000, ])

ppc_intervals(y = stan_data_NJ$nonzero_positives,
              yrep = y_rep[1:1000, ])

fit_summary_NJ <- summary(fit_model_NJ)

medians_rt <- fit_summary_NJ$summary[, '50%'][230 :(230+226) ]
min_rt_50_interval <- fit_summary_NJ$summary[, '25%'][230 :(230+226) ]
max_rt_50_interval <- fit_summary_NJ$summary[, '75%'][230 :(230+226) ]
min_rt_95_interval <- fit_summary_NJ$summary[, '2.5%'][230 :(230+226) ]
max_rt_95_interval <- fit_summary_NJ$summary[, '97.5%'][230 :(230+226) ]


ggplot(data = NULL, aes(x = NJ_data$date, y = medians_rt)) + 
  geom_line() + 
  xlab('Date') +
  ylab('') +
  ggtitle( 'NJ r_t')+
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  geom_vline(xintercept = NJ_data$date[1]) +
  geom_ribbon(aes(ymin = min_rt_50_interval, ymax = max_rt_50_interval), alpha= 0.5, fill = 'darkred') +
  geom_ribbon(aes(ymin = min_rt_95_interval, ymax = max_rt_95_interval), alpha= 0.1, fill = 'darkred')


