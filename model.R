library(rstan)
source('data_cleaning_us.R')
source('utility_functions.R')

### model for New Jersey

X <- get_model_data(data, 'NJ', Sys.Date() - 2)

p_delay <- get_delay_distribution()

nonzero_days <- which(X$total != 0)

stan_data <- list(N = nrow(X),
                  conv_gt = get_gt_convolution(nrow(X)),
                  length_delay = length(p_delay),
                  p_delay = p_delay,
                  exposures = exposures_from_total(X$total),
                  N_nonzero = length(nonzero_days),
                  nonzero_positives = X$positive[nonzero_days],
                  nonzero_days = nonzero_days
  
)


compiled_model <- stan_model('rt_model.stan')

fit_model <- sampling(compiled_model, data=stan_data, iter = 1000)





