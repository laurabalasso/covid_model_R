library(ramify)
source('../R/regional_model_data.R')


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


## Dummies for red/orange/yellow area:
## generation of a matrix for each level of risk (yellow = medium , orange = medium-high , red = high)
## given the selected regions and dates

yellow_dummies <- function(Y_df, regions, dates){
  Y <- Y_df
  Y <- Y %>% filter(date %in% dates)
  Y_0 <- matrix(0, nrow = length(dates[dates<as.Date('2020-11-06')]), ncol = length(regions))
  Y <- rbind(Y_0, as.matrix(Y[regions]))
  return(Y)
}

orange_dummies <- function(O_df, regions, dates){
  O <- O_df
  O <- O %>% filter(date %in% dates)
  O_0 <- matrix(0, nrow = length(dates[dates<as.Date('2020-11-06')]), ncol = length(regions))
  O <- rbind(O_0, as.matrix(O[regions]))
  return(O)
}

red_dummies <- function(R_df, regions, dates){
  R <- R_df
  R <- R %>% filter(date %in% dates)
  R_0 <- matrix(0, nrow = length(dates[dates<as.Date('2020-11-06')]), ncol = length(regions))
  R <- rbind(R_0, as.matrix(R[regions]))
  
  return(R)
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

get_hier_data <- function(hier_data, regions, Y, O, R, initial_date =  min(hier_data$date) , final_date = Sys.Date() - 2 ,buffer_days = 10){
  positives <- get_hier_positives(hier_data,regions,  initial_date, final_date, buffer_days)
  total <- get_hier_totals(hier_data,regions, initial_date, final_date, buffer_days)
  exposures <- get_hier_exposures(total)
  dates <- seq(initial_date - buffer_days, final_date, by = 'day')
  nonzero_days <- which(apply(total, 1, prod) != 0)
  nonzero_dates <- dates[nonzero_days]
  yellow_dummies <- yellow_dummies(Y,regions, dates)
  orange_dummies <- orange_dummies(O,regions, dates)
  red_dummies <- red_dummies(R,regions, dates)
  
  l <- list(positives = positives,total = total, exposures = exposures, dates = dates, 
            nonzero_days = nonzero_days, nonzero_dates = nonzero_dates, yellow_dummies = yellow_dummies,
            orange_dummies = orange_dummies, red_dummies=red_dummies)
  
  return(l)
  
}


