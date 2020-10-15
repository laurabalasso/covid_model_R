## delay distribution

get_delay_distribution <- function(){
  p_delay <- read.csv('./data/p_delay.csv')
  p_delay <- p_delay$p_delay
}

## generation time distribution

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


get_model_data <- function(data, region, date, start_date =  data$date[1] ,buffer_days = 10){
  df <- data[data$region == region & data$date <= as.Date(date) & data$date >= as.Date(start_date) , ]
  df <- df[-1]
  
  first_idx <- which(df$positive != 0)[1]
  df <- df[first_idx : nrow(df), ] 
  
  buffer_idx <- seq(df$date[1] - buffer_days , df$date[1] -1, by = "day")
  buffer <- data.frame(date = buffer_idx, positive =  rep(0, buffer_days),total = rep(0,buffer_days))
  df <- rbind(buffer, df, make.row.names = FALSE)
  
  
  return(df)
}



### compute convolution matrix from generation time distribution

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



