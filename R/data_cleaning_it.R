library(dplyr)
library(magrittr)
library(plyr)
library(tidyr)
library(ramify)

#############################
###   total = tamponi     ###
#############################

## data_it <- read.csv('../data/data_regions.csv')

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


negatives <- data.frame(region = data_it$region[which(data_it$positive<0)], date =data_it$date[which(data_it$positive<0)])


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




### From 2020-11-06 regions grouped in 3 areas based on risk level 
## (yellow = medium , orange = medium-high , red = high)

d1 <- as.Date('2020-11-06')

r_d1 <- c('Calabria', 'Lombardia', 'Piemonte', 'Valle d\'Aosta')
o_d1 <-c('Puglia', 'Sicilia')
y_d1 <- c('Abruzzo', 'Basilicata', 'Campania', 'Emilia-Romagna', "Friuli Venezia Giulia",
                  'Lazio', 'Liguria','Marche', 'Molise', 'P.A. Bolzano', 'P.A. Trento', 'Sardegna',
                  'Toscana', 'Umbria', 'Veneto')


d2 <- as.Date('2020-11-11')
r_d2 <- c('P.A. Bolzano')
o_d2 <- c('Abruzzo', 'Basilicata', 'Liguria' , 'Toscana', 'Umbria' )
y_d2 <- c()

d3 <- as.Date('2020-11-15')
r_d3<-c('Campania', 'Toscana')
o_d3 <- c('Emilia-Romagna', 'Friuli Venezia Giulia', 'Marche')
y_d3 <- c()

d4 <- as.Date('2020-11-22')
r_d4 <- c('Abruzzo')
o_d4 <- c()
y_d4 <- c()


d5 <- as.Date('2020-11-29')
y_d5 <- c('Liguria', 'Sicilia')
o_d5 <- c('Calabria', 'Lombardia' , 'Piemonte')
r_d5 <- c()

d6 <- as.Date('2020-12-06')
y_d6 <- c('Emilia-Romagna', 'Friuli Venezia Giulia', 'Marche', 'Puglia', 'Umbria' )
o_d6 <- c('Campania', 'Toscana', 'Valle D\'Aosta' , 'P.A. Bolzano')
r_d6 <- c()

ordinanze <- c(d1, d2, d3, d4, d5, d6)
in_yellow <- list(y_d1, y_d2, y_d3, y_d4, y_d5, y_d6)
in_orange <- list(o_d1, o_d2, o_d3, o_d4, o_d5, o_d6)
in_red <-    list(r_d1, r_d2, r_d3, r_d4, r_d5, r_d6)


## generation of a matrix for each level of risk 
## each column refers to a region
## takes value 1 if the region has that level of risk on that date (0 otherwise)

Y<-O<-R <- data.frame(date=unique(data_it$date[data_it$date >= '2020-11-06']))

for(r in unique(data_it$region)){
  Y[r] <- O[r] <- R[r] <- rep(0, nrow(Y))
}


for(r in 1:length(unique(data_it$region))){
  reg <- unique(data_it$region)[r]
  for(i in 1:length(ordinanze)){
    if(reg %in% unlist(in_yellow[i]))
      Y[which(Y$date >= ordinanze[i]), r+1] <- rep(1, length(which(Y$date >= ordinanze[i])))
    else if(reg %in% unlist(in_orange[i]) | reg %in% unlist(in_red[i]))
      Y[which(Y$date >= ordinanze[i]), r+1] <- rep(0, length(which(Y$date >= ordinanze[i])))
  }
}


for(r in 1:length(unique(data_it$region))){
  reg <- unique(data_it$region)[r]
  for(i in 1:length(ordinanze)){
    if(reg %in% unlist(in_orange[i]))
      O[which(O$date >= ordinanze[i]), r+1] <- rep(1, length(which(O$date >= ordinanze[i])))
    else if(reg %in% unlist(in_yellow[i]) | reg %in% unlist(in_red[i]))
      O[which(O$date >= ordinanze[i]), r+1] <- rep(0, length(which(O$date >= ordinanze[i])))
  }
}

for(r in 1:length(unique(data_it$region))){
  reg <- unique(data_it$region)[r]
  for(i in 1:length(ordinanze)){
    if(reg %in% unlist(in_red[i]))
      R[which(R$date >= ordinanze[i]), r+1] <- rep(1, length(which(R$date >= ordinanze[i])))
    else if(reg %in% unlist(in_yellow[i]) | reg %in% unlist(in_orange[i]))
      R[which(R$date >= ordinanze[i]), r+1] <- rep(0, length(which(R$date >= ordinanze[i])))
  }
}


