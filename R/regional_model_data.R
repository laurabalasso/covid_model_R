library(ramify)
source('../R/auxiliary_distributions.R')


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

get_gt_convolution_ln1 <- function(N_days){
  gt <- get_generation_time_distribution_ln1()
  convolution_ready_gt = matrix(rep(0, N_days*(N_days-1)), N_days - 1, N_days )
  
  for(t in 1:(N_days-1)){
    begin <- max(1, t - length(gt) + 2)
    slice_update <- rev(gt[2 : (t - begin + 2)]) 
    convolution_ready_gt[t , begin : (begin + length(slice_update) - 1)] <- slice_update 
  }
  
  return(convolution_ready_gt)
}


get_gt_convolution_ln2 <- function(N_days){
  gt <- get_generation_time_distribution_ln2()
  convolution_ready_gt = matrix(rep(0, N_days*(N_days-1)), N_days - 1, N_days )
  
  for(t in 1:(N_days-1)){
    begin <- max(1, t - length(gt) + 2)
    slice_update <- rev(gt[2 : (t - begin + 2)]) 
    convolution_ready_gt[t , begin : (begin + length(slice_update) - 1)] <- slice_update 
  }
  
  return(convolution_ready_gt)
}

get_gt_convolution_g <- function(N_days){
  gt <- get_generation_time_distribution_g()
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



