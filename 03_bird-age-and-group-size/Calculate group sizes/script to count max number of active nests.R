
################################################################################
# Script to calculate number of active breeding pairs during first broods
# using Google Sheet nest check data
# Heather Kenny-Duddela
# June 3, 2024
################################################################################

# libraries
library(tidyverse)
library(ggplot2)

# load data
mb <- read.csv("Make Believe 2022_nest check.csv")

# get rid of first and third rows which just indicates which dates were checked
mb2 <- mb[-c(1,3),]

# check cells for "CI" and "X"
active.index <- as.data.frame(which(mb2 == "CI" | 
                                      mb2 == "X", arr.ind = T))

mb.active <- mb2[c(1,unique(active.index$row)),]
colnames(mb.active)[2:155] <- mb.active[1,2:155]

# for each active nest, pull out CI dates and X dates

# storage
dates <- data.frame(nest=NA, CI.date=NA, X.date=NA)

# loop to look through each row

nests <- length(mb.active$Make)-1

for (i in 2:nests) {
  # get CI dates
  index.ci <- which(mb.active[i,] == "CI")
  dates.ci <- colnames(mb.active)[index.ci]
  # get inactive dates
  index.x <- which(mb.active[i,] == "X" | mb.active[i,]=="x")
  dates.x <- colnames(mb.active[index.x])
  # number of needed rows
  size <- length(index.ci)
  
  # fill in storage
  dates.temp <- data.frame(nest=rep(NA,size), CI.date=rep(NA,size), 
                           X.date=rep(NA,size))
  dates.temp$nest <- mb.active[i,1]
  dates.temp$CI.date <- dates.ci
  dates.temp$X.date <- dates.x
  
  # combine with external storage
  dates <- rbind(dates, dates.temp)
}

## Mark active periods for each nest

# dates April 16 to Sep 16 (154 days)

library(lubridate)

days <- ymd("2022-09-16") - ymd("2022-04-16")

# remove extra row from dates
dates2 <- dates[-1,]
dates2$CI.date <- as.Date(as.factor(dates2$CI.date), format="%d-%b-%y")
dates2$X.date <- as.Date(as.factor(dates2$X.date), format="%d-%b-%y")

# make matrix where columns are days and rows are nests
nest.mat <- as.data.frame(matrix(nrow=nests, ncol=155, 0))
colnames(nest.mat) <- c("nest", as.character(seq(ymd("2022-04-16"),
                                                 ymd("2022-09-16"),by=1)))
nest.mat[,1] <- mb.active$Make[2:11]

# fill in cells with 1s for active times for each nest

for (i in 1:length(dates2$nest)) {
  # get start and end dates
  period <- dates2[i, ]
  # get vector of column names (dates)
  mat.colnames <- colnames(nest.mat)
  # get column numbers
  start <- which(mat.colnames == as.character(period$CI.date))
  stop <- which(mat.colnames == as.character(period$X.date))
  # fill in row of nest.mat
  nest.mat[nest.mat$nest == period$nest, start:stop] <- 1
  
}

# get column sums
ColSums <- colSums(nest.mat[,2:155])
# add to nest matrix
nest.mat.sum <- rbind(nest.mat, c(NA, ColSums))

toplot <- data.frame(date = ymd(colnames(nest.mat[-1])),
                        active = ColSums)

# bar plot to see when max number of active nests occur
ggplot(toplot, aes(x=date, y=active)) + geom_point() + geom_line() +
  geom_vline(xintercept = ymd("2022-06-20"), color="red") +
  ggtitle("Make Believe active nests across time 2022")

ggsave("Make Believe active nests over time 2022.png", h=4, w=5)

mb.max <- max(ColSums) # 7 in 1st brood and 7 max

#-------------------------------------------------------------------------------
# Try same process with CHR

chr <- read.csv("Colorado Horse Rescue 2022_nest checks.csv")

# get rid of second row which just indicates which dates were checked
chr2 <- chr[-c(2),]

# check cells for "CI" and "X"
chr.active.index <- as.data.frame(which(chr2 == "CI" | 
                                      chr2 == "X" |
                                        chr2=="x", arr.ind = T))

chr.active <- chr2[c(1,unique(chr.active.index$row)),]
colnames(chr.active)[2:155] <- chr.active[1,2:155]

# remove problem nests, but include them in final count
# nest 36 became active on 7/13/2022
# but then couldn't be checked for the rest of the season
chr.active2 <- subset(chr.active, chr.active$Checked.!="36")


# for each active nest, pull out CI dates and X dates
# storage
chr.dates <- data.frame(nest=NA, CI.date=NA, X.date=NA)

# loop to look through each row

nests.chr <- length(chr.active2$Checked.)-1

for (i in 2:nests) {
  # get CI dates
  index.ci <- which(chr.active2[i,] == "CI")
  dates.ci <- colnames(chr.active2)[index.ci]
  # get inactive dates
  index.x <- which(chr.active2[i,] == "X" | chr.active2[i,]=="x")
  dates.x <- colnames(chr.active2[index.x])
  # number of needed rows
  size <- length(index.ci)
  
  # fill in storage
  dates.temp <- data.frame(nest=rep(NA,size), CI.date=rep(NA,size), 
                           X.date=rep(NA,size))
  dates.temp$nest <- chr.active2[i,1]
  dates.temp$CI.date <- dates.ci
  dates.temp$X.date <- dates.x
  
  # combine with external storage
  chr.dates <- rbind(chr.dates, dates.temp)
}


## Mark active periods for each nest

# dates April 16 to Sep 16 (154 days, inclusive)

# remove extra row from dates
chr.dates2 <- chr.dates[-1,]
chr.dates2$CI.date <- as.Date(as.factor(chr.dates2$CI.date), format="%d-%b-%y")
chr.dates2$X.date <- as.Date(as.factor(chr.dates2$X.date), format="%d-%b-%y")

# make matrix where columns are days and rows are nests
chr.nest.mat <- as.data.frame(matrix(nrow=nests.chr, ncol=155, 0))
colnames(chr.nest.mat) <- c("nest", as.character(seq(ymd("2022-04-16"),
                                                 ymd("2022-09-16"),by=1)))
chr.nest.mat[,1] <- chr.active2$Checked.[2:60]

# fill in cells with 1s for active times for each nest

for (i in 1:length(chr.dates2$nest)) {
  # get start and end dates
  period <- chr.dates2[i, ]
  # get vector of column names (dates)
  mat.colnames <- colnames(chr.nest.mat)
  # get column numbers
  start <- which(mat.colnames == as.character(period$CI.date))
  stop <- which(mat.colnames == as.character(period$X.date))
  # fill in row of nest.mat
  chr.nest.mat[chr.nest.mat$nest == period$nest, start:stop] <- 1
}

# get column sums
chr.ColSums <- colSums(chr.nest.mat[,2:155])
# add to nest matrix
chr.nest.mat.sum <- rbind(chr.nest.mat, c(NA, chr.ColSums))

chr.toplot <- data.frame(date = ymd(colnames(chr.nest.mat[-1])),
                     active = chr.ColSums)

# bar plot to see when max number of active nests occur
ggplot(chr.toplot, aes(x=date, y=active)) + geom_point() + geom_line() +
  geom_vline(xintercept = ymd("2022-06-20"), color="red") +
  ggtitle("CHR active nests over time 2022")

ggsave("CHR active nests over time 2022.png", h=4, w=5)

chr.max <- max(chr.ColSums) # 33 active nests in 1st brood max, 33 max in 2nd
# broods including problem nest 


#-------------------------------------------------------------------------------
# Calculate number of max nests for Blue Cloud

# load Blue Cloud data
bc <- read.csv("Blue Cloud 2022_nest check.csv")

# delete extra row two
bc2 <- bc[-c(2),]

# check cells for "CI" and "X"
bc.active.index <- as.data.frame(which(bc2 == "CI" | 
                                          bc2 == "X" |
                                          bc2=="x", arr.ind = T))

bc.active <- bc2[c(1,unique(bc.active.index$row)),]
colnames(bc.active)[2:154] <- bc.active[1,2:154]

# for each active nest, pull out CI dates and X dates
# storage
bc.dates <- data.frame(nest=NA, CI.date=NA, X.date=NA)

# loop to look through each row

nests.bc <- length(bc.active$Checked.)-1

for (i in 2:nests.bc) {
  # get CI dates
  index.ci <- which(bc.active[i,] == "CI")
  dates.ci <- colnames(bc.active)[index.ci]
  # get inactive dates
  index.x <- which(bc.active[i,] == "X" | bc.active[i,]=="x")
  dates.x <- colnames(bc.active[index.x])
  # number of needed rows
  size <- length(index.ci)
  
  # fill in storage
  dates.temp <- data.frame(nest=rep(NA,size), CI.date=rep(NA,size), 
                           X.date=rep(NA,size))
  dates.temp$nest <- bc.active[i,1]
  dates.temp$CI.date <- dates.ci
  dates.temp$X.date <- dates.x
  
  # combine with external storage
  bc.dates <- rbind(bc.dates, dates.temp)
}


## Mark active periods for each nest

# dates April 16 to Sep 15 (153 days, inclusive)

# remove extra row from dates
bc.dates2 <- bc.dates[-1,]
bc.dates2$CI.date <- as.Date(as.factor(bc.dates2$CI.date), format="%d-%b-%y")
bc.dates2$X.date <- as.Date(as.factor(bc.dates2$X.date), format="%d-%b-%y")

# make matrix where columns are days and rows are nests
bc.nest.mat <- as.data.frame(matrix(nrow=nests.bc, ncol=154, 0))
colnames(bc.nest.mat) <- c("nest", as.character(seq(ymd("2022-04-16"),
                                                     ymd("2022-09-15"),by=1)))
bc.nest.mat[,1] <- bc.active$Checked.[2:13]

# fill in cells with 1s for active times for each nest

for (i in 1:length(bc.dates2$nest)) {
  # get start and end dates
  period <- bc.dates2[i, ]
  # get vector of column names (dates)
  mat.colnames <- colnames(bc.nest.mat)
  # get column numbers
  start <- which(mat.colnames == as.character(period$CI.date))
  stop <- which(mat.colnames == as.character(period$X.date))
  # fill in row of nest.mat
  bc.nest.mat[bc.nest.mat$nest == period$nest, start:stop] <- 1
}

# get column sums
bc.ColSums <- colSums(bc.nest.mat[,2:154])
# add to nest matrix
bc.nest.mat.sum <- rbind(bc.nest.mat, c(NA, bc.ColSums))

bc.toplot <- data.frame(date = ymd(colnames(bc.nest.mat[-1])),
                         active = bc.ColSums)

# bar plot to see when max number of active nests occur
ggplot(bc.toplot, aes(x=date, y=active)) + geom_point() + geom_line() +
  geom_vline(xintercept = ymd("2022-06-20"), color="red") +
  ggtitle("Blue Cloud active nests over time 2022")

ggsave("Blud Cloud active nests over time 2022.png", h=4, w=5)

bc.max <- max(bc.ColSums) # BC 6 active nests in 1st broods, 7 active nest max
# 7 pairs total - 3 in main barn, one in back stall, one in trailer, one in gazebo,
# one in wash room

#-------------------------------------------------------------------------------
# Calculate max number of active nests for Cooks

# load cooks data
cook <- read.csv("Cooks 2022_nest check.csv")

# delete extra row two
cook2 <- cook[-c(2),]

# check cells for "CI" and "X"
cook.active.index <- as.data.frame(which(cook2 == "CI" | 
                                         cook2 == "X" |
                                         cook2=="x", arr.ind = T))

cook.active <- cook2[c(1,unique(cook.active.index$row)),]
colnames(cook.active)[2:155] <- cook.active[1,2:155]

# for each active nest, pull out CI dates and X dates
# storage
cook.dates <- data.frame(nest=NA, CI.date=NA, X.date=NA)

# loop to look through each row

nests.cook <- length(cook.active$Checked.)-1

for (i in 2:nests.cook) {
  # get CI dates
  index.ci <- which(cook.active[i,] == "CI")
  dates.ci <- colnames(cook.active)[index.ci]
  # get inactive dates
  index.x <- which(cook.active[i,] == "X" | cook.active[i,]=="x")
  dates.x <- colnames(cook.active[index.x])
  # number of needed rows
  size <- length(index.ci)
  
  # fill in storage
  dates.temp <- data.frame(nest=rep(NA,size), CI.date=rep(NA,size), 
                           X.date=rep(NA,size))
  dates.temp$nest <- cook.active[i,1]
  dates.temp$CI.date <- dates.ci
  dates.temp$X.date <- dates.x
  
  # combine with external storage
  cook.dates <- rbind(cook.dates, dates.temp)
}


## Mark active periods for each nest

# dates April 16 to Sep 22 (160 days, inclusive)

# remove extra row from dates
cook.dates2 <- cook.dates[-1,]
cook.dates2$CI.date <- as.Date(as.factor(cook.dates2$CI.date), format="%d-%b-%y")
cook.dates2$X.date <- as.Date(as.factor(cook.dates2$X.date), format="%d-%b-%y")

# make matrix where columns are days and rows are nests
cook.nest.mat <- as.data.frame(matrix(nrow=nests.cook, ncol=161, 0))
colnames(cook.nest.mat) <- c("nest", as.character(seq(ymd("2022-04-16"),
                                                    ymd("2022-09-22"),by=1)))
cook.nest.mat[,1] <- cook.active$Checked.[2:16]

# fill in cells with 1s for active times for each nest

for (i in 1:length(cook.dates2$nest)) {
  # get start and end dates
  period <- cook.dates2[i, ]
  # get vector of column names (dates)
  mat.colnames <- colnames(cook.nest.mat)
  # get column numbers
  start <- which(mat.colnames == as.character(period$CI.date))
  stop <- which(mat.colnames == as.character(period$X.date))
  # fill in row of nest.mat
  cook.nest.mat[cook.nest.mat$nest == period$nest, start:stop] <- 1
}

# get column sums
cook.ColSums <- colSums(cook.nest.mat[,2:161])
# add to nest matrix
cook.nest.mat.sum <- rbind(cook.nest.mat, c(NA, cook.ColSums))

cook.toplot <- data.frame(date = ymd(colnames(cook.nest.mat[-1])),
                        active = cook.ColSums)

# bar plot to see when max number of active nests occur
ggplot(cook.toplot, aes(x=date, y=active)) + geom_point() + geom_line() +
  geom_vline(xintercept = ymd("2022-06-20"), color="red") +
  ggtitle("Cooks active nests over time 2022") + ylim(0,11)

ggsave("Cooks active nests over time 2022.png", h=4, w=5)

cook.max <- max(cook.ColSums) # 7 active nests in 1st brood, 9 total

#-------------------------------------------------------------------------------
# Calculate group size for Struthers

# load data
st <- read.csv("Struthers 2022_nest check.csv")

# delete extra row two
st2 <- st[-c(2),]

# check cells for "CI" and "X"
st.active.index <- as.data.frame(which(st2 == "CI" | 
                                           st2 == "X" |
                                           st2=="x", arr.ind = T))

st.active <- st2[c(1,unique(st.active.index$row)),]
colnames(st.active)[2:155] <- st.active[1,2:155]

# for each active nest, pull out CI dates and X dates
# storage
st.dates <- data.frame(nest=NA, CI.date=NA, X.date=NA)

# loop to look through each row

nests.st <- length(st.active$Checked.)-1

for (i in 2:nests.st) {
  # get CI dates
  index.ci <- which(st.active[i,] == "CI")
  dates.ci <- colnames(st.active)[index.ci]
  # get inactive dates
  index.x <- which(st.active[i,] == "X" | st.active[i,]=="x")
  dates.x <- colnames(st.active[index.x])
  # number of needed rows
  size <- length(index.ci)
  
  # fill in storage
  dates.temp <- data.frame(nest=rep(NA,size), CI.date=rep(NA,size), 
                           X.date=rep(NA,size))
  dates.temp$nest <- st.active[i,1]
  dates.temp$CI.date <- dates.ci
  dates.temp$X.date <- dates.x
  
  # combine with external storage
  st.dates <- rbind(st.dates, dates.temp)
}


## Mark active periods for each nest

# dates April 16 to Sep 16 (154 days, inclusive)

# remove extra row from dates
st.dates2 <- st.dates[-1,]
st.dates2$CI.date <- as.Date(as.factor(st.dates2$CI.date), format="%d-%b-%y")
st.dates2$X.date <- as.Date(as.factor(st.dates2$X.date), format="%d-%b-%y")

# make matrix where columns are days and rows are nests
st.nest.mat <- as.data.frame(matrix(nrow=nests.st, ncol=155, 0))
colnames(st.nest.mat) <- c("nest", as.character(seq(ymd("2022-04-16"),
                                                      ymd("2022-09-16"),by=1)))
st.nest.mat[,1] <- st.active$Checked.[2:6]

# fill in cells with 1s for active times for each nest

for (i in 1:length(st.dates2$nest)) {
  # get start and end dates
  period <- st.dates2[i, ]
  # get vector of column names (dates)
  mat.colnames <- colnames(st.nest.mat)
  # get column numbers
  start <- which(mat.colnames == as.character(period$CI.date))
  stop <- which(mat.colnames == as.character(period$X.date))
  # fill in row of nest.mat
  st.nest.mat[st.nest.mat$nest == period$nest, start:stop] <- 1
}

# get column sums
st.ColSums <- colSums(st.nest.mat[,2:155])
# add to nest matrix
st.nest.mat.sum <- rbind(st.nest.mat, c(NA, st.ColSums))

st.toplot <- data.frame(date = ymd(colnames(st.nest.mat[-1])),
                          active = st.ColSums)

# bar plot to see when max number of active nests occur
ggplot(st.toplot, aes(x=date, y=active)) + geom_point() + geom_line() +
  geom_vline(xintercept = ymd("2022-06-20"), color="red") +
  ggtitle("Struthers active nests over time 2022")

ggsave("Struthers active nests over time 2022.png", h=4, w=5)

st.max <- max(st.ColSums) # 2 in 1st broods, 3 total

#-------------------------------------------------------------------------------
# Group size for Urban Farm Girlz

# Only 2 active nests, 5 birds banded, 2 pairs. No need to look at overlap of active nests


