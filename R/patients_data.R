

download.file("https://github.com/beoutbreakprepared/nCoV2019/raw/master/latest_data/latestdata.tar.gz", 'data/patients_data.tar.gz')
untar('data/patients_data.tar.gz', exdir = 'data')

patients <- read.csv(pipe("cut -f6,10,12 -d',' data/latestdata.csv"), header = TRUE)


# There's an errant reversed date

patients$date_onset_symptoms[which(patients$date_onset_symptoms =="01.31.2020") ] <- "31.01.2020"


## keep only observations with both symptom onset and test confirmation date 

patients <- patients[ - which(patients$date_onset_symptoms == ''), ]
patients <- patients[ - which(patients$date_confirmation == ''), ]


## check valid date and convert to standard date format

IsDate <- function(mydate, date.format = "%d.%m.%y") {
  tryCatch(!is.na(as.Date(mydate, date.format)),  
           error = function(err) {FALSE})  
}


patients <- patients[which(IsDate(patients$date_onset_symptoms) & IsDate(patients$date_confirmation)), ]

patients$date_onset_symptoms <- as.Date(patients$date_onset_symptoms, format = '%d.%m.%y')
patients$date_confirmation <- as.Date(patients$date_confirmation , format ='%d.%m.%y' )


# Only keep records where confirmed > onset

patients <- patients[which(patients$date_onset_symptoms < patients$date_confirmation), ]

# Mexico has many cases that are all confirmed on the same day regardless
# of onset date, so we filter it out.

patients <- patients[which(patients$country != 'Mexico'), ]
  

### delays 

max_delay<- 60

delays <- patients$date_confirmation - patients$date_onset_symptoms

delays <- delays[delays < max_delay]

## get distribution

incubation_days = 5

p_delay <- table(delays)

p_delay <- p_delay / sum(p_delay)

delay <-data.frame(p_delay = c( rep(0, incubation_days), as.vector(p_delay)))

write.csv(delay, 'data/p_delay.csv')

  