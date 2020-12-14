library(rstan)
library(bayesplot)
source('R/data_cleaning_it.R')
source('R/hierarchical_model_data.R')
source('R/rt_plots.R')

## regions <- unique(data_it$region)

#regions <- c('Lazio', 'Lombardia', 'Abruzzo', 'Veneto', 'Molise', 'Basilicata')

regions <- c('Lombardia', 'Lazio')
hier_data <- get_hier_data(data_it, regions, initial_date = as.Date('2020-09-15'))


### lockdown effect

lockdown_start <- as.Date('2020-03-09')
lockdown_end <- as.Date('2020-05-03')

lockdown <- rep(0, length(hier_data$dates))
lockdown[which(hier_data$dates> lockdown_start+10 & hier_data$dates <= lockdown_end)] <- 1
grow <- which(hier_data$dates>=lockdown_start &hier_data$dates <= lockdown_start +10)
lockdown[grow] <- (grow - which(hier_data$dates ==lockdown_start))^2 /100


p_delay <- get_delay_distribution()

stan_data_hier <- list(J = length(regions),
                       N = nrow(hier_data$exposures),
                       N_nonzero = length(hier_data$nonzero_days),
                       nonzero_days = hier_data$nonzero_days,
                       conv_gt = get_gt_convolution_ln2(nrow(hier_data$exposures)),
                       length_delay = length(p_delay),
                       p_delay = p_delay,
                       exposures = hier_data$exposures,
                       nonzero_positives = hier_data$positives[hier_data$nonzero_days ,]
)


compiled_hier <- stan_model('stan/hier_rt_model.stan')

fit_hier <- sampling(compiled_hier, data = stan_data_hier, iter= 100, cores= getOption("mc.cores", 1L))

print(fit_hier, pars = 'y_rep')


mcmc_trace(as.array(fit_hier, pars = c('tau')),
           np = nuts_params(fit_hier)
)
mcmc_trace(as.array(fit_hier, pars = c('phi')),
           np = nuts_params(fit_hier)
)


mcmc_intervals(fit_hier, pars=c( "mu[1]","mu[2]", "mu[3]", "mu[4]", "mu[5]", "mu[6]", "mu[7]",
                                     "mu[8]","mu[9]", "mu[10]", "mu[11]", "mu[12]", "mu[13]", "mu[14]",
                                     "mu[15]","mu[16]", "mu[17]", "mu[18]", "mu[19]", 
                                 "mu[20]","mu[21]", "mu[22]", "mu[23]", "mu[24]",  "mu[25]","mu[26]", 
                                 "mu[27]","mu[28]", "mu[29]", "mu[30]", "mu[31]", "mu[32]", "mu[33]",
                                 "mu[34]","mu[35]", "mu[36]", "mu[37]", "mu[38]", "mu[39]", "mu[40]",
                                 "mu[41]","mu[42]", "mu[43]", "mu[44]", "mu[45]","mu[46]","mu[48]", 
                                 "mu[49]","mu[50]", "mu[51]", "mu[52]", "mu[53]", "mu[54]", "mu[55]",
                                 "mu[56]", "mu[57]", "mu[58]", "mu[59]"))

mcmc_hist(fit_hier, pars = c('tau'), binwidth = 0.01)

### Effective sample size 

ratios1 <- neff_ratio(fit_hier, pars = c('mu'))
mcmc_neff(ratios1) + yaxis_text(hjust = 1)

ratios2 <- neff_ratio(fit_hier, pars = c('r_t'))
mcmc_neff(ratios2) + yaxis_text(hjust = 1)

check_n_eff(fit_hier)

parcoord_with_divs <- mcmc_parcoord(
  as.array(fit_hier, pars = c("tau", "mu[1]", "mu[2]", "mu[3]", "mu[4]")),
  np = nuts_params(fit_hier)
)
parcoord_with_divs

scatter_with_divs <- mcmc_scatter(
  as.array(fit_hier),
  pars = c("mu[4]", 'tau'),
  transform = list('tau' = "log"),
  np = nuts_params(fit_hier)
)
scatter_with_divs


regional_yrep_idx <- function(region, regions_vector, nonzero_days){
  region_idx <- which(regions_vector == region)
  yrep_idx <- (region_idx-1)* length(nonzero_days) + 1
  range <- yrep_idx : (yrep_idx + length(nonzero_days)-1)
  return(range)
}



y_rep <- as.matrix(fit_hier, pars = "y_rep")
ppc_dens_overlay(y = as.vector(stan_data_hier$nonzero_positives), y_rep[1:1000, ]) 
ppc_dens_overlay(y = stan_data_hier$nonzero_positives[, which(regions == 'Lazio')], y_rep[1:1000,regional_yrep_idx('Lazio', regions, stan_data_hier$nonzero_days)]) 


ppc_intervals(
  y = stan_data_hier$nonzero_positives[, which(regions == 'Lazio')],
  yrep = y_rep[, regional_yrep_idx('Lazio', regions,stan_data_hier$nonzero_days)]
)



#### Grouped plots


groups <- function(regions, nonzero_days){
  group <- rep(regions[1], length(nonzero_days))
  for(r in 2:length(regions))
    group <- c(group, rep(regions[r], length(nonzero_days)))
  
  return(group)
}


ppc_stat(y=as.vector(stan_data_hier$nonzero_positives), yrep =y_rep,stat="mean")
ppc_stat_grouped(y=as.vector(stan_data_hier$nonzero_positives), yrep =y_rep, group = groups(regions, stan_data_hier$nonzero_days) ,stat="mean", binwidth=0.5)


mean_y_rep<-colMeans(y_rep)
std_resid<-(as.vector(stan_data_hier$nonzero_positives)-mean_y_rep)/sqrt(mean_y_rep)
qplot(mean_y_rep,std_resid)+hline_at(2)+hline_at(-2)


ppc_intervals_grouped(y=as.vector(stan_data_hier$nonzero_positives),yrep=y_rep,group=groups(regions, stan_data_hier$nonzero_days))+labs(x="date",y="new positive cases")


### R_t curve


plot_rt_hier(hier_data, fit_hier, regions, 'Lombardia')

