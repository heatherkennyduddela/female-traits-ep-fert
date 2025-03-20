################################################################################
# Script to calculate estimated age of birds from 2022 using longitudinal database
# Heather Kenny-Duddela
# May 29, 2024
################################################################################

### libraries
library(dplyr)
library(ggplot2)


### load data

# 2022 adult banding data
band <- read.csv("2022 Adult Data_2022-12-27.csv")

# adult database
adult <- read.csv("AdultData2008-2023_withColor.csv")
# nestling database
nestling <- read.csv("NestlingData2008-2022.csv")

# pull out instances of 2022 adults from adult database
# since database includes banding through 2023, all 2022 adults will show up
adult2 <- adult[adult$band %in% band$band, ]

# pull out instances of 2022 adults from nestling banding
nestling2 <- nestling[nestling$band %in% band$band, ]

# only keep band, date, site, sex
adult3 <- adult2 %>% select(band, site, date_captured, sex, recap)
colnames(adult3)[3] <- "date"

nestling3 <- nestling2 %>% select(band, site, date)
nestling3$sex <- NA
nestling3$recap <- NA

# combine nestling and adult resights
both <- rbind(adult3, nestling3)

# for each bird, get first sighting and calculate estimated age
age <- both %>% group_by(band) %>%
  summarise(first_date = min(date), 
            first_year = as.numeric(substr(first_date, 1, 4)),
            index = which(date == first_date),
            adult.nestl = sex[index])

# check for duplicated bands
dup <- age[duplicated(age$band), ]
# remove duplicate
age2 <- age[-which(duplicated(age$band)), ]

# calculate age as number of years breeding. So birds first seen in 2022 get 1
# instead of zero
age$est.age <- 2023 - age$first_year
# for birds first caught as nestlings, they are zero in the first year, so 
# subtract one from the estimated age. Unknown sex is usual a fledgling
age$est.age[which(is.na(age$adult.nestl) | age$adult.nestl=="unknown")] <- 
  age$est.age[which(is.na(age$adult.nestl)|
                      age$adult.nestl=="unknown")] -1

# add site from 2022
age2 <- left_join(age, select(band, band, site), by="band")
colnames(age2)[7] <- "site_2022"

# add column for age certainty score
# "exact" means banded as neslting, "good" means the site previously had a high
# banding effort for adults, "low" means the site was not fully banded in previous
# years so unbanded birds could be older than 1

# check adult banding from 2021
prev.year <- subset(adult, adult$year==2021)
# count up adults from each site in 2020
prev.year.sum <- prev.year %>% group_by(site) %>%
  summarise(adults.band = n())

age2$certainty <- "low"
age2$certainty[which(age2$site_2022 == "Make Believe" |
                       age2$site_2022 == "Grizz" |
                       age2$site_2022 == "Blue Cloud" |
                       age2$site_2022 == "Schaap " |
                       age2$site_2022 == "Mary Ann's" |
                       age2$site_2022 == "Struthers" |
                       age2$site_2022 == "Cooks" |
                       age2$site_2022 == "Urban Farm Girlz" |
                       age2$site_2022 == "Cathy's" )] <- "good"
age2$certainty[which(age2$first_year<=2021)] <- "good" # known ASY
age2$certainty[which(is.na(age2$adult.nestl) |
                             age2$adult.nestl=="unknown")] <- "exact"


# histogram of ages colored by certainty
ggplot(age2, aes(x=est.age, fill=certainty)) + 
  geom_histogram(binwidth = 1, color="black") +
  scale_x_continuous(name="Estimated age (years)", limits=c(0,10), breaks=c(1:10)) +
  ggtitle("Histogram of estimated ages for 2022 adults\ncolored by certainty")

ggsave("adult ages 2022 histogram.png")

# save table of ages
write.csv(age2, "table of estimated ages 2022.csv")

