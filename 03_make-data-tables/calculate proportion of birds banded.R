
# Script to calculate proportion of birds banded based on site sizes
# Heather Kenny-Duddela
# March 13, 2025

library(dplyr)
library(lubridate)
library(tidyr)

# load data, original input file
band22 <- read.csv("input-files/AdultData2022_fromDatabase.csv")

#-------------------------------------------------------------------------------

# check for fledglings
unique(band22$sex)
fledge <- filter(band22, sex=="unknown" | is.na(sex))

adult <- filter(band22, sex!="unknown" & !is.na(sex))
# remove ? from sex
adult$sex2 <- substr(adult$sex, 1, 1)

# summarize females and males per site
adult.summary <- adult %>% group_by(site, sex2) %>%
  summarise(n=n())

# keep only core sites
adult.core <- filter(adult.summary, site=="Blue Cloud" | site=="Cathys" |
                       site=="Colorado Horse Rescue" | site=="Cooks" |
                       site=="Dome House" | site=="Karens" |
                       site=="Make Believe" | site=="McCauley" |
                       site=="Struthers" | site=="Urban Farm Girlz")


# convert to wide
adult.wide <- pivot_wider(adult.core, names_from = sex2, values_from = n)
# correct Urban Farm Girlz
adult.wide[10,2] <- 2

adult.wide$total <- adult.wide$F + adult.wide$M

# add site sizes (alphabetical)
site.size <- c(7, 1, 33, 9, 1, 1, 7, 2, 3, 2)

adult.wide$est.site.size <- site.size*2

# adjust site size based on number of banded birds and known UNB birds
# CHR had at least 1 UNB female, so total females=38, and there were two
# polygamous males, so 38 + 36 = 74 for adj.size. It also had at least 3 UNB
# males, which makes sense give the number of banded males. 

# Cooks had at one UNB female, so adj.size=21

# Blue Cloud had one UNB male, so keep est.size=14

# Make Believe had zero UNB, but an extra male arrived later in the season. 
# Keep adj.size=15
adult.wide$adj.size <- adult.wide$est.site.size
adult.wide$adj.size[which(adult.wide$site=="Colorado Horse Rescue")] <- 74
adult.wide$adj.size[which(adult.wide$site=="Cooks")] <- 21
adult.wide$adj.size[which(adult.wide$site=="Make Believe")] <- 15

# add column for proportion banded
adult.wide$prop.band <- adult.wide$total / adult.wide$adj.size

# overall prop banded: 
sum(adult.wide$total)/sum(adult.wide$adj.size) # 0.9513889, 137 out of 144


