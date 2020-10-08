## delay distribution


get_delay_distribution <- function(){
  p_delay <- read.csv('data/p_delay.csv')
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



