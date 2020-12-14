library(rstan)
library(bayesplot)
source('R/data_cleaning_it.R')
source('R/regional_model_data.R')
source('R/rt_plots.R')

### Model for Lombardia with school opening slope

data_lombardia <- get_model_data(data_it, 'Lombardia', initial_date = as.Date('2020-03-01'), final_date = as.Date('2020-05-06'))

## school effect
ggplot(data_lombardia, aes(x = date, y = positive)) + geom_line() + 
  geom_vline( xintercept = as.Date('2020-09-14') + 10, color = 'red')

school_opening <- as.Date('2020-09-14')

### gradual increasing school effect 
school<-c()

for(i in 1:length(data_lombardia$date)) {
  if (data_lombardia$date[i]<=school_opening) { school[i]=0}
  else if(data_lombardia$date[i]>=school_opening+10) {school[i]=1}
  else  {school[i]=( (i - sum(data_lombardia$date <=school_opening ))^2)/100}
}

### lockdown effect

lockdown_start <- as.Date('2020-03-09')
lockdown_end <- as.Date('2020-05-03')

lockdown <- rep(0, length(data_lombardia$date))
lockdown[which(data_lombardia$date> lockdown_start+10 & data_lombardia$date <= lockdown_end)] <- 1
grow <- which(data_lombardia$date>=lockdown_start & data_lombardia$date <= lockdown_start +10)
lockdown[grow] <- (grow - which(data_lombardia$date ==lockdown_start))^2 /100
  


## dummy

school <- as.numeric(data_lombardia$date > (school_opening + 10))

p_delay <- get_delay_distribution()

nonzero_days_l <- which(data_lombardia$total != 0)

stan_data_school <- list(N = nrow(data_lombardia),
                    conv_gt = get_gt_convolution_ln2(nrow(data_lombardia)),
                    length_delay = length(p_delay),
                    p_delay = p_delay,
                    exposures = exposures_from_total(data_lombardia$total),
                    N_nonzero = length(nonzero_days_l),
                    nonzero_positives = data_lombardia$positive[nonzero_days_l],
                    nonzero_days = nonzero_days_l,
                    lockdown = lockdown
)


compiled_model_school <- stan_model('rt_model_schools.stan')

fit_model_lomb_school <- sampling(compiled_model_school, data=stan_data_school, iter = 2000)


print(fit_model_lomb_school, pars = 'r_t')
print(fit_model_lomb_school, pars = 'beta_school')


mcmc_hist(fit_model_lomb_school, pars = 'beta_school')


## trace plots

mcmc_trace(as.array(fit_model_lomb_school, pars = c('r_t[1]', 'r_t[10]', 'r_t[50]', 'r_t[100]')),
           np = nuts_params(fit_model_lomb_school)
)

mcmc_trace(as.array(fit_model_lomb_school, pars = c('mu[1]', 'mu[10]', 'mu[50]', 'mu[100]')),
           np = nuts_params(fit_model_lomb_school)
)


mcmc_trace(
  as.array(fit_model_lomb_school, pars = c('alpha', 'seed')),
  np = nuts_params(fit_model_lomb_school)
)




## effective sample size ratios 

ratios1 <- neff_ratio(fit_model_lomb_school, pars = c('mu'))
mcmc_neff(ratios1) + yaxis_text(hjust = 1)

ratios2 <- neff_ratio(fit_model_lomb_school, pars = c('r_t'))
mcmc_neff(ratios2) + yaxis_text(hjust = 1)

check_n_eff(fit_model_lomb_school)



## posterior predictive check

y_rep <- as.matrix(fit_model_lomb_school, pars = "y_rep")
ppc_dens_overlay(y = data_lombardia$positive[nonzero_days_l], y_rep[1:1000, ]) 

ppc_intervals(
  y = data_lombardia$positive[nonzero_days_l],
  yrep = y_rep
)


ppc_stat(y = stan_data_school$nonzero_positives, yrep = y_rep, stat = 'mean')


### R_t curve

plot_rt(data_lombardia, fit_model_lomb_school, 'Lombardia')