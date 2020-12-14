## delay distribution 

get_delay_distribution <- function(){
  p_delay <- read.csv('../data/p_delay.csv')
  p_delay <- p_delay$p_delay
}

## generation time distribution 
# Source: https://www.ijidonline.com/article/S1201-9712(20)30119-3/pdf 

get_generation_time_distribution_ln1 <- function(){
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

# Source https://www.frontiersin.org/articles/10.3389/fphy.2020.00347/full 
get_generation_time_distribution_ln2 <- function(){
  mean_si = 4.9
  std_si = 4.4
  mu_si = log(mean_si ** 2 / sqrt(std_si ** 2 + mean_si ** 2))
  sigma_si = sqrt(log(std_si ** 2 / mean_si ** 2 + 1))
  
  x <- 0:19
  y <- plnorm(x, mu_si, sigma_si)
  y <- c(0, diff(y))
  y <- y/sum(y)
  
  return(y)
}

get_generation_time_distribution_g <- function(){
  shape <-  1.87
  rate <- 0.28
  
  x <- 0:19
  y <- pgamma(x, shape = shape, rate = rate)
  y <- c(0, diff(y))
  y <- y/sum(y)
  
  return(y)
}


