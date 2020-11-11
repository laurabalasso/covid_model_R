library(ramify)

## delay distribution 

get_delay_distribution <- function(){
  p_delay <- read.csv('data/p_delay.csv')
  p_delay <- p_delay$p_delay
}

## generation time distribution 
# Source: https://www.ijidonline.com/article/S1201-9712(20)30119-3/pdf 

get_generation_time_distribution <- function(){
  mean_si = 4.7
  std_si = 2.9
  mu_si = log(mean_si ** 2 / sqrt(std_si ** 2 + mean_si ** 2))
  sigma_si = sqrt(log(std_si ** 2 / mean_si ** 2 + 1))
  
  x <- 0:19
  y <- plnorm(x, mu_si, sigma_si)
  y <- c(0, diff(y))
  y <- y/sum(y)
  
  return(y)
}


##################################################
### Preprocessing data for single region model ###
##################################################

## Generate dataset for single region model fitting
# with columns date, positive (daily new positives) and totals (daily swabs increment)

## Arguments:
# data = complete national data containg daily values of new positives and swabs increment for each region
# region  = name of the region/state to select (type = string)
# initial_date / final_date = boundaries for the time period to model (type = Date) 
# buffer_days =  number of blank days to pad on the leading edge of the time series 

get_model_data <- function(data, region,initial_date =  data$date[1], final_date = Sys.Date() -2 ,buffer_days = 10){
  df <- data[data$region == region & data$date <= as.Date(final_date) & data$date >= as.Date(initial_date) , ]
  df <- df[-1]
  
  first_idx <- which(df$positive != 0)[1]
  df <- df[first_idx : nrow(df), ] 
  
  buffer_idx <- seq(df$date[1] - buffer_days , df$date[1] -1, by = "day")
  buffer <- data.frame(date = buffer_idx, positive =  rep(0, buffer_days),total = rep(0,buffer_days))
  df <- rbind(buffer, df, make.row.names = FALSE)
  
  
  return(df)
}



### compute an auxiliary convolution matrix from generation time distribution (for fitting speedup)

get_gt_convolution <- function(N_days){
  gt <- get_generation_time_distribution()
  convolution_ready_gt = matrix(rep(0, N_days*(N_days-1)), N_days - 1, N_days )
  
  for(t in 1:(N_days-1)){
    begin <- max(1, t - length(gt) + 2)
    slice_update <- rev(gt[2 : (t - begin + 2)]) 
    convolution_ready_gt[t , begin : (begin + length(slice_update) - 1)] <- slice_update 
  }
  
  return(convolution_ready_gt)
}



### compute exposures from totals

exposures_from_total <- function(total){
  exposures = clip(total, max(total)*0.1, 1e09)
}


###############################################
## Preprocessing data for hierarchical model ##
###############################################

## New positives matrix
## each column is a daily time series of new positives for a single region 

get_hier_positives <- function(hier_data, regions, initial_date =  min(hier_data$date) , final_date = Sys.Date() - 2 ,buffer_days = 10){
  M <- matrix(NA, nrow = as.numeric(final_date - initial_date) + buffer_days +1 , ncol = length(regions))
  for(j in 1:length(regions)){
    r <- get_model_data(hier_data, regions[j], initial_date, final_date, buffer_days )
    
    diff <-  nrow(M) - length(r$positive) 
    if(diff > 0)
      r <- c(rep(0, diff ), r$positive)
    else
      r <- r$positive
    
    M[, j] <- r
  }
  return(M)
}


## Total swabs daily increment matrix
## each column is a daily time series of swabs increment for a single region 

get_hier_totals <- function(hier_data, regions, initial_date =  min(hier_data$date) , final_date = Sys.Date() - 2 ,buffer_days = 10){
  M <- matrix(NA, nrow = as.numeric(final_date - initial_date) + buffer_days +1 , ncol = length(regions))
  for(j in 1:length(regions)){
    r <- get_model_data(hier_data, regions[j], initial_date, final_date, buffer_days )
    
    diff <-  nrow(M) - length(r$total) 
    if(diff > 0)
      r <- c(rep(0, diff ), r$total)
    else
      r <- r$total
    
    M[, j] <- r
  }
  return(M)
}


## Exposures matrix from total swabs time series 

get_hier_exposures <- function(total){
  M <- matrix(NA, nrow = nrow(total) , ncol = ncol(total))
  for(j in 1:ncol(M)){
    M[, j] <- clip(total[, j], max(total[, j])*0.1, 1e09)
  }
  return(M)
}


###  Generate necessary data for hierarchical model fitting ##

## Arguments: 
# hier_data = complete national data containg daily values of new positives and swabs increment for each region
# regions = vector of regions' names to include in the model (type = string)
# initial_date / final_date = boundaries for the time period to model (type = Date) 
# buffer_days =  number of blank days to pad on the leading edge of the time series 

## Returns a list containg:
# positives: matrix of daily new positives for each selected region (cols = regions) 
# total: matrix of daily increment of swabs for each selected region (cols = regions) 
# exposures: matrix of exposures for each selected region (cols = regions) for test volume correction of new positives 
# dates: vector of selected dates
# nonzero_days: indices for which swabs increment is > 0 for all selected regions 
# nonzero_dates: dates for which swabs increment is > 0 for all selected regions 

# If the first available date is not the same for all regions, some of them will have a longer buffer

get_hier_data <- function(hier_data, regions, initial_date =  min(hier_data$date) , final_date = Sys.Date() - 2 ,buffer_days = 10){
  positives <- get_hier_positives(hier_data,regions,  initial_date, final_date, buffer_days)
  total <- get_hier_totals(hier_data,regions, initial_date, final_date, buffer_days)
  exposures <- get_hier_exposures(total)
  dates <- seq(initial_date - buffer_days, final_date, by = 'day')
  nonzero_days <- which(apply(total, 1, prod) != 0)
  nonzero_dates <- dates[nonzero_days]
  
  l <- list(positives = positives,total = total, exposures = exposures, dates = dates, 
            nonzero_days = nonzero_days, nonzero_dates = nonzero_dates)
  
  return(l)
  
}


######################
### model checking ###
######################

## Function to check for problematic effective sample size values
# source: https://github.com/betanalpha/knitr_case_studies/blob/b71d65d44731ce90bbfc769f3cbc8355efac262f/divergences_and_bias/stan_utility.R#L47

check_n_eff <- function(fit) {
  fit_summary <- summary(fit, probs = c(0.5))$summary
  N <- dim(fit_summary)[[1]]
  
  iter <- dim(rstan::extract(fit)[[1]])[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    ratio <- fit_summary[,5][n] / iter
    if (ratio < 0.001) {
      print(sprintf('n_eff / iter for parameter %s is %s!',
                    rownames(fit_summary)[n], ratio))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('n_eff / iter looks reasonable for all parameters')
  else
    print('  n_eff / iter below 0.001 indicates that the effective sample size has likely been overestimated')
}




