library(ggplot2)

#####################
##### Rt plots  #####
#####################


### Retrieve rt medians and intervals for a given model and region

rt_intervals <- function(data, fit){
  fit_summary <- summary(fit)
  
  rt_idx <- which(rownames(fit_summary$summary) == 'r_t[1]')
  rt_median <- fit_summary$summary[rt_idx: (rt_idx + nrow(data) - 1), '50%']
  min_rt_50_interval <- fit_summary$summary[rt_idx: (rt_idx + nrow(data) - 1), '25%']
  max_rt_50_interval <- fit_summary$summary[rt_idx: (rt_idx + nrow(data) - 1), '75%']
  min_rt_95_interval <- fit_summary$summary[rt_idx: (rt_idx + nrow(data) - 1), '2.5%']
  max_rt_95_interval <- fit_summary$summary[rt_idx: (rt_idx + nrow(data) - 1), '97.5%']
  
  return(list(median = rt_median, min_50 = min_rt_50_interval, max_50 = max_rt_50_interval, min_95 = min_rt_95_interval, max_95 = max_rt_95_interval))
  
}

rt_intervals_hier <- function(data, fit, regions_vector, region){
  summary <- summary(fit)
  region_idx = which(regions_vector == region)
  rt_1 <- which(rownames(summary$summary) == paste('r_t[1,', region_idx,  ']', sep=''))
  rt_index <- seq(rt_1, rt_1 + length(data$dates) * length(regions_vector)  -1 , by=length(regions_vector))
  rt_median <- summary$summary[rt_index, '50%']
  min_rt_50_interval <- summary$summary[rt_index, '25%']
  max_rt_50_interval <- summary$summary[rt_index, '75%']
  min_rt_95_interval <- summary$summary[rt_index, '2.5%']
  max_rt_95_interval <- summary$summary[rt_index, '97.5%']
  
  return(list(median = rt_median, min_50 = min_rt_50_interval, max_50 = max_rt_50_interval, min_95 = min_rt_95_interval, max_95 = max_rt_95_interval))
}

### Plot rt curve for a single region

plot_rt <- function(data, fit, region){
  intervals <- rt_intervals(data, fit)
  p<- ggplot(data = NULL, aes(x = data$date, y = intervals$median)) + 
    geom_line() + 
    xlab('Date') +
    ylab('') +
    ggtitle( paste(region, ' r_t'))+
    geom_hline(yintercept=1, linetype="dashed", color = "red") +
    geom_vline(xintercept = data$date[1]) +
    geom_ribbon(aes(ymin = intervals$min_50, ymax = intervals$max_50), alpha= 0.5, fill = 'darkred') +
    geom_ribbon(aes(ymin = intervals$min_95, ymax = intervals$max_95), alpha= 0.1, fill = 'darkred')
  
  return(p)
  
}

plot_rt_hier <- function(data, fit, regions_vector, region){
  
  intervals <- rt_intervals_hier(data, fit, regions_vector, region)
  
  p <- ggplot(  ) + 
    geom_line(aes(x = data$dates, y=intervals$median ), color='darkred')+
    xlab('Date') +
    ylab('') +
    geom_hline(yintercept=1, linetype="dashed", color = "red") +
    geom_vline(xintercept = data$dates[1]) +
    geom_ribbon(aes(x = data$dates, ymin = intervals$min_50, ymax = intervals$max_50), alpha= 0.5, fill = 'darkred') +
    geom_ribbon(aes(x = data$dates,ymin =intervals$min_95, ymax = intervals$max_95), alpha= 0.1, fill = 'darkred')
  
  return(p)
}



get_plot_limits <- function(plot) {
  gb = ggplot_build(plot)
  xmin = gb$layout$panel_params[[1]]$x.range[1]
  xmax = gb$layout$panel_params[[1]]$x.range[2]
  ymin = gb$layout$panel_params[[1]]$y.range[1]
  ymax = gb$layout$panel_params[[1]]$y.range[2]
  list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}

### overlay with rt plot of the single region model 
plot_overlay <- function(plot,data_2, fit_2, regions,  region){
  
  intervals_2 <- rt_intervals_hier(data_2, fit_2, regions, region)
  
  p <- plot + geom_line(aes(x=data_2$dates, y = intervals_2$median), color = 'navyblue') +
    geom_ribbon(aes(x = data_2$dates,ymin = intervals_2$min_50, ymax = intervals_2$max_50), alpha= 0.4, fill = 'royalblue1') +
    geom_ribbon(aes(x = data_2$dates,ymin = intervals_2$min_95, ymax = intervals_2$max_95), alpha= 0.1, fill = 'royalblue1')+
    ggtitle(paste('Rt', region , sep = ' '))
  
  return(p)
  
}
