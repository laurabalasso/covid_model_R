library(rstan)
library(bayesplot)
source('R/data_cleaning_it.R')
source('R/regional_model_data.R')
source('R/rt_plots.R')

data_fvg <- get_model_data(data_it, 'Friuli Venezia Giulia', initial_date = as.Date('2020-08-30'), final_date = Sys.Date()-2)

p_delay <- get_delay_distribution()

nonzero_days_fvg <- which(data_fvg$total != 0)

stan_data_fvg <- list(N = nrow(data_fvg),
                    conv_gt = get_gt_convolution_ln2(nrow(data_fvg)),
                    length_delay = length(p_delay),
                    p_delay = p_delay,
                    exposures = exposures_from_total(data_fvg$total),
                    N_nonzero = length(nonzero_days_fvg),
                    nonzero_positives = data_fvg$positive[nonzero_days_fvg],
                    nonzero_days = nonzero_days_fvg
)


compiled_model_fvg <- stan_model('stan/rt_model.stan')

fit_model_fvg <- sampling(compiled_model_fvg, data=stan_data_fvg, iter = 100, cores=getOption("mc.cores", 1L))


y_rep <- as.matrix(fit_model_fvg, pars = "y_rep")
ppc_dens_overlay(y = data_fvg$positive[nonzero_days_fvg], y_rep[1:100, ]) 

ppc_intervals(
  y = data_fvg$positive[nonzero_days_fvg],
  yrep = y_rep
)



ppc_stat(y = stan_data_fvg$nonzero_positives, yrep = y_rep, stat = 'mean')


### R_t curve

plot_rt(data_fvg, fit_model_fvg, 'FVG')


