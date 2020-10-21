library(dplyr)
library(magrittr)
library(plyr)
library(tidyr)
library(ramify)

#############################
###   total = tamponi     ###
#############################


data_it <- read.csv(url("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni.csv"))

data_it$data <- as.Date(data_it$data)


data_it <- arrange(data_it, denominazione_regione, data)

data_it <- subset(data_it, select = c(denominazione_regione, data, nuovi_positivi, tamponi))


## Transform cumulative total tests to daily increase

data_it$tamponi <- c(NaN, diff(data_it$tamponi))

data_it <- data_it %>% drop_na() 

data_it$tamponi <- clip(data_it$tamponi, 0 ,NULL)

colnames(data_it)[colnames(data_it) == 'denominazione_regione'] <- 'region'
colnames(data_it)[colnames(data_it) == 'data'] <- 'date'
colnames(data_it)[colnames(data_it) == 'nuovi_positivi'] <- 'positive'
colnames(data_it)[colnames(data_it) == 'tamponi'] <- 'total'


#############################
### total = casi testati  ###
#############################

data_it_ct <- read.csv(url("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni.csv"))

data_it_ct$data <- as.Date(data_it_ct$data)


data_it_ct <- arrange(data_it_ct, denominazione_regione, data)

data_it_ct <- subset(data_it_ct, select = c(denominazione_regione, data, nuovi_positivi, casi_testati))


## Transform cumulative total tests to daily increase

data_it_ct$casi_testati <- c(NaN, diff(data_it_ct$casi_testati))

data_it_ct <- data_it_ct %>% drop_na() 

data_it_ct$casi_testati <- clip(data_it_ct$casi_testati, 0 ,NULL)

colnames(data_it_ct)[colnames(data_it_ct) == 'denominazione_regione'] <- 'region'
colnames(data_it_ct)[colnames(data_it_ct) == 'data'] <- 'date'
colnames(data_it_ct)[colnames(data_it_ct) == 'nuovi_positivi'] <- 'positive'
colnames(data_it_ct)[colnames(data_it_ct) == 'casi_testati'] <- 'total'



