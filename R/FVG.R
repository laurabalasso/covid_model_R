library(rstan)
library(bayesplot)
library(ggplot2)
source('R/data_cleaning_it.R')
source('R/utility_functions.R')

data_fvg <- get_model_data(data_it, 'Friuli Venezia Giulia', initial_date = as.Date('2020-07-30'), final_date = Sys.Date())

p_delay <- get_delay_distribution()

nonzero_days_fvg <- which(data_fvg$total != 0)

stan_data_fvg <- list(N = nrow(data_fvg),
                    conv_gt = get_gt_convolution(nrow(data_fvg)),
                    length_delay = length(p_delay),
                    p_delay = p_delay,
                    exposures = exposures_from_total(data_fvg$total),
                    N_nonzero = length(nonzero_days_fvg),
                    nonzero_positives = data_fvg$positive[nonzero_days_fvg],
                    nonzero_days = nonzero_days_fvg
)


compiled_model_fvg <- stan_model('stan/rt_model.stan')

fit_model_fvg <- sampling(compiled_model_fvg, data=stan_data_fvg, iter = 3000, cores=getOption("mc.cores", 1L))


y_rep <- as.matrix(fit_model_fvg, pars = "y_rep")
ppc_dens_overlay(y = data_fvg$positive[nonzero_days_fvg], y_rep[1:1000, ]) 

ppc_intervals(
  y = data_fvg$positive[nonzero_days_fvg],
  yrep = y_rep
)



ppc_stat(y = stan_data_fvg$nonzero_positives, yrep = y_rep, stat = 'mean')


### R_t curve

fit_summary_fvg<- summary(fit_model_fvg)

rt_idx <- which(rownames(fit_summary_fvg$summary) == 'r_t[44]')
medians_rt <- fit_summary_fvg$summary[rt_idx: (rt_idx + stan_data_fvg$N - 44), '50%']
min_rt_50_interval <- fit_summary_fvg$summary[rt_idx: (rt_idx + stan_data_fvg$N - 44), '25%']
max_rt_50_interval <- fit_summary_fvg$summary[rt_idx: (rt_idx + stan_data_fvg$N - 44), '75%']
min_rt_95_interval <- fit_summary_fvg$summary[rt_idx: (rt_idx + stan_data_fvg$N - 44), '2.5%']
max_rt_95_interval <- fit_summary_fvg$summary[rt_idx: (rt_idx + stan_data_fvg$N -44), '97.5%']


ggplot(data = NULL, aes(x = data_fvg$date[44:115], y = medians_rt)) + 
  geom_line() + 
  xlab('Date') +
  ylab('') +
  ggtitle( 'Rt FVG') +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  geom_vline(xintercept = data_fvg$date[44]) +
  geom_ribbon(aes(ymin = min_rt_50_interval, ymax = max_rt_50_interval), alpha= 0.5, fill = 'darkred') +
  geom_ribbon(aes(ymin = min_rt_95_interval, ymax = max_rt_95_interval), alpha= 0.1, fill = 'darkred')






