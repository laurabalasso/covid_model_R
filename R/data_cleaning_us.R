library(dplyr)
library(magrittr)
library(plyr)
library(tidyr)
library(ramify) ## clip function


data<- read.csv(url("https://api.covidtracking.com/v1/states/daily.csv"))

colnames(data)[colnames(data) == 'state'] <- 'region'

data$date <-as.Date(as.character(data$date),format="%Y%m%d")
data <- arrange(data, region , date )

data <- subset(data, select = c(region, date, positive, total))


## pochi dati per alcune regioni 
data <- data[which(data$region != 'GU'), ]
data <- data[which(data$region != 'MP'), ]
data <- data[which(data$region != 'AS'), ]
data <- data[which(data$region != 'PR'), ]
data <- data[which(data$region != 'VI'), ]

# On Jun 5 Covidtracking started counting probable cases too
# which increases the amount by 5014.
# https://covidtracking.com/screenshots/MI/MI-20200605-184320.png

data$positive[data$region == 'MI' & data$date >= as.Date('2020-06-05')] <- data$positive[data$region == 'MI' & data$date >= as.Date('2020-06-05')] - 5014


# From CT: On June 19th, LDH removed 1666 duplicate and non resident cases
# after implementing a new de-duplicaton process.
data$positive[data$region == 'LA' & data$date >= as.Date('2020-06-19')] <- data$positive[data$region == 'LA' & data$date >= as.Date('2020-06-19')] + 1666
data$total[data$region == 'LA' & data$date >= as.Date('2020-06-19')] <- data$total[data$region == 'LA' & data$date >= as.Date('2020-06-19')] + 1666


## transform cumulative positives and totals to daily increase
data$positive <- c(NaN, diff(data$positive))
data$total <- c(NaN, diff(data$total))

data <- data %>% drop_na() 

data$positive <- clip(data$positive, 0, NULL)
data$total <- clip(data$total, 0 ,NULL)


# Michigan missed 6/18 totals and lumped them into 6/19 so we've
# divided the totals in two and equally distributed to both days.

data$total[data$region == 'MI' & data$date == as.Date("2020-06-18")] <- 1471
data$total[data$region == 'MI' & data$date == as.Date("2020-06-19")] <- 1471



# Note that when we set total to zero, the model ignores that date. See
# the likelihood function in GenerativeModel.build

data$positive[data$region == 'NJ' & data$date == as.Date("2020-05-11") ] <- 0 
data$total[data$region == 'NJ'& data$date == as.Date("2020-05-11")] <- 0


data$positive[data$region == 'NJ' & data$date == as.Date("2020-07-25") ] <- 0 
data$total[data$region == 'NJ'& data$date == as.Date("2020-07-25")] <- 0


data$positive[data$region == 'CA' & data$date == as.Date("2020-04-22") ] <- 0 
data$total[data$region == 'CA'& data$date == as.Date("2020-04-22")] <- 0

data$positive[data$region == 'SC' & data$date == as.Date("2020-06-26") ] <- 0 
data$total[data$region == 'SC'& data$date == as.Date("2020-06-26")] <- 0

data$positive[data$region == 'OR' & data$date == as.Date("2020-06-28") ] <- 174
data$total[data$region == 'OR'& data$date == as.Date("2020-06-28")] <- 3296

data$positive[data$region == 'OH' & data$date == as.Date("2020-07-01") ] <- 0 
data$total[data$region == 'OH'& data$date == as.Date("2020-07-01")] <- 0
data$positive[data$region == 'OH' & data$date == as.Date("2020-07-09") ] <- 0 
data$total[data$region == 'OH'& data$date == as.Date("2020-07-09")] <- 0

data$positive[data$region == 'NV' & data$date == as.Date("2020-07-02") ] <- 0 
data$total[data$region == 'NV'& data$date == as.Date("2020-07-02")] <- 0

data$positive[data$region == 'AL' & data$date == as.Date("2020-07-09") ] <- 0 
data$total[data$region == 'AL'& data$date == as.Date("2020-07-09")] <- 0

data$positive[data$region == 'AR' & data$date == as.Date("2020-07-10") ] <- 0 
data$total[data$region == 'AR'& data$date == as.Date("2020-07-10")] <- 0

data$positive[data$region == 'MS' & data$date == as.Date("2020-07-12") ] <- 0 
data$total[data$region == 'MS'& data$date == as.Date("2020-07-12")] <- 0

data$positive[data$region == 'CT' & data$date == as.Date("2020-07-17") ] <- 0 
data$total[data$region == 'CT'& data$date == as.Date("2020-07-17")] <- 0
data$positive[data$region == 'CT' & data$date == as.Date("2020-07-21") ] <- 0 
data$total[data$region == 'CT'& data$date == as.Date("2020-07-21")] <- 0

data$positive[data$region == 'DC' & data$date == as.Date("2020-08-04") ] <- 0 
data$total[data$region == 'DC'& data$date == as.Date("2020-08-04")] <- 0


data$positive[data$region == 'PA' & data$date == as.Date("2020-06-03") ] <- 0 
data$total[data$region == 'PA'& data$date == as.Date("2020-06-03")] <- 0
data$positive[data$region == 'PA' & data$date == as.Date("2020-04-21") ] <- 0 
data$total[data$region == 'PA'& data$date == as.Date("2020-04-21")] <- 0
data$positive[data$region == 'PA' & data$date == as.Date("2020-05-20") ] <- 0 
data$total[data$region == 'PA'& data$date == as.Date("2020-05-20")] <- 0

data$positive[data$region == 'HI' & data$date == as.Date("2020-08-07") ] <- 0 
data$total[data$region == 'HI'& data$date == as.Date("2020-08-07")] <- 0

data$positive[data$region == 'TX' & data$date == as.Date("2020-08-08") ] <- 0 
data$total[data$region == 'TX'& data$date == as.Date("2020-08-08")] <- 0
data$positive[data$region == 'TX' & data$date == as.Date("2020-08-11") ] <- 0 
data$total[data$region == 'TX'& data$date == as.Date("2020-08-11")] <- 0

data$positive[data$region == 'DE' & data$date == as.Date("2020-08-14") ] <- 0 
data$total[data$region == 'DE'& data$date == as.Date("2020-08-14")] <- 0

data$positive[data$region == 'SD' & data$date == as.Date("2020-08-26") ] <- 0 
data$total[data$region == 'SD'& data$date == as.Date("2020-08-26")] <- 0

data$positive[data$region == 'WA' & data$date >= as.Date("2020-09-22") & data$date <= as.Date("2020-09-24") ] <- 0 
data$total[data$region == 'WA'& data$date >= as.Date("2020-09-22") & data$date <= as.Date("2020-09-24")] <- 0


# Zero out any rows where positive tests equal or exceed total reported tests
# Do not act on Wyoming as they report positive==total most days

filtering_date <- as.Date('2020-07-27')
zero_filter <- (data$positive >= data$total) & (data$date >= filtering_date) & (data$region != 'WY')


data$positive[zero_filter] <- 0
data$total[zero_filter] <- 0


