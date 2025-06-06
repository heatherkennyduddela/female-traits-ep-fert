
################################################################################
# Correlation tables for raw data, within sex, between mates, with RS
# Heather Kenny-Duddela
# Oct 28, 2024
################################################################################

# libraries
library(lubridate)
library(dplyr)
library(Hmisc) # for correlation tables

# load data
fem <- read.csv("output-files/female_season_table_2025-03-13.csv")
male <- read.csv("output-files/male_season_table_2025-03-13.csv")

# add julian dates
fem$ci_1_julian <- yday(fem$ci_1)
male$socF_ci_julian <- yday(male$socF_ci)

fem$ci_1 <- ymd(fem$ci_1)
male$socF_ci <- ymd(male$socF_ci)

# change factor order
fem$age.class <- factor(fem$age.class, levels=c("SY","ASY"))
male$age.class <- factor(male$age.class, levels=c("SY","ASY"))



#-------------------------------------------------------------------------------
# Correlations within females
#-------------------------------------------------------------------------------

## Include wing, tail, color, ci, total RS, fecundity

# take numeric columns only
fem.cor <- fem %>%
  select(mean_rwl, tail, throat_avg_bright, breast_avg_bright,
                  belly_avg_bright, vent_avg_bright, ci_1_julian, tot.chick,
                  tot_eggs)

# calculate correlations
fem.table <- rcorr(as.matrix(fem.cor), type="spearman")
# pull out correlation values
fem.R <- fem.table$r
# pull out p-vals
fem.p <- fem.table$P

## Define notions for significance levels; spacing is important.
mystars <- ifelse(fem.p < .0001, "****", 
                  ifelse(fem.p < .001, "*** ",
                         ifelse(fem.p < .01, "**  ",
                         ifelse(fem.p < .05, "*   ", "    "))))

## truncate the correlation matrix to two decimal
fem.R2 <- format(round(cbind(rep(-1.11, ncol(fem.cor)), fem.R), 2))[,-1]

## build a new matrix that includes the correlations with their appropriate stars
Rnew <- matrix(paste(fem.R2, mystars, sep=""), ncol=ncol(fem.cor))
diag(Rnew) <- paste(diag(fem.R), " ", sep="")
rownames(Rnew) <- colnames(fem.cor)
colnames(Rnew) <- paste(colnames(fem.cor), "", sep="")

# Hide upper triangle
fem.upper<-Rnew
fem.upper[upper.tri(Rnew)]<-""
fem.upper<-as.data.frame(fem.upper)

write.csv(fem.upper, "output-files/female corr table_2025-03-13.csv")

#-------------------------------------------------------------------------------
# Correlations within males
#-------------------------------------------------------------------------------

# take numeric columns only
male.cor <- select(male, mean_rwl, tail, throat_avg_bright, breast_avg_bright,
                  belly_avg_bright, vent_avg_bright, socF_ci_julian, tot.chick,
                  season.prop.wp)

# calculate correlations
male.table <- rcorr(as.matrix(male.cor), type="spearman")
# pull out correlation values
male.R <- male.table$r
# pull out p-vals
male.p <- male.table$P

## Define notions for significance levels; spacing is important.
mystars <- ifelse(male.p < .0001, "****", 
                  ifelse(male.p < .001, "*** ",
                         ifelse(male.p < .01, "**  ",
                         ifelse(male.p < .05, "*   ", "    "))))

## truncate the correlation matrix to two decimal
male.R2 <- format(round(cbind(rep(-1.11, ncol(male.cor)), male.R), 2))[,-1]

## build a new matrix that includes the correlations with their appropriate stars
Rnew <- matrix(paste(male.R2, mystars, sep=""), ncol=ncol(male.cor))
diag(Rnew) <- paste(diag(male.R), " ", sep="")
rownames(Rnew) <- colnames(male.cor)
colnames(Rnew) <- paste(colnames(male.cor), "", sep="")

# Hide upper triangle
male.upper<-Rnew
male.upper[upper.tri(Rnew)]<-""
male.upper<-as.data.frame(male.upper)

write.csv(male.upper, "output-files/male corr table_2025-03-13.csv")

#-------------------------------------------------------------------------------
# correlations between social pairs - from female table
#-------------------------------------------------------------------------------

## Include wing, tail, color, ci, total RS

# take numeric columns only
fem.socM.cor <- select(fem, mean_rwl, tail, throat_avg_bright, breast_avg_bright,
                  belly_avg_bright, vent_avg_bright, ci_1_julian, socM_rwl,
                  socM_tail, socM_t.avg.bright, socM_r.avg.bright,
                  socM_b.avg.bright, socM_v.avg.bright, tot.chick)

# calculate correlations
fem.socM.table <- rcorr(as.matrix(fem.socM.cor), type="spearman")
# pull out correlation values
fem.socM.R <- fem.socM.table$r
# pull out p-vals
fem.socM.p <- fem.socM.table$P

## Define notions for significance levels; spacing is important.
mystars <- ifelse(fem.socM.p < .0001, "****", 
                  ifelse(fem.socM.p < .001, "*** ",
                         ifelse(fem.socM.p < .01, "**  ",
                         ifelse(fem.socM.p < .05, "*   ", "    "))))

## truncate the correlation matrix to two decimal
fem.socM.R2 <- format(round(cbind(rep(-1.11, ncol(fem.socM.cor)), fem.socM.R), 2))[,-1]

## build a new matrix that includes the correlations with their appropriate stars
Rnew <- matrix(paste(fem.socM.R2, mystars, sep=""), ncol=ncol(fem.socM.cor))
diag(Rnew) <- paste(diag(fem.socM.R), " ", sep="")
rownames(Rnew) <- colnames(fem.socM.cor)
colnames(Rnew) <- paste(colnames(fem.socM.cor), "", sep="")

# Hide upper triangle
fem.socM.upper<-Rnew
fem.socM.upper[upper.tri(Rnew)]<-""
fem.socM.upper<-as.data.frame(fem.socM.upper)

write.csv(fem.socM.upper, "output-files/fem socM corr table_2025-03-13.csv")

#-------------------------------------------------------------------------------
# correlations between social pairs - from male table
#-------------------------------------------------------------------------------

# take numeric columns only
male.socF.cor <- select(male, socF_rwl, socF_tail, socF_t.avg.bright, 
                   socF_r.avg.bright, socF_b.avg.bright, socF_v.avg.bright,
                   socF_ci_julian, mean_rwl, tail, throat_avg_bright, 
                   breast_avg_bright,
                   belly_avg_bright, vent_avg_bright,  tot.chick,
                   season.prop.wp)

# calculate correlations
male.socF.table <- rcorr(as.matrix(male.socF.cor), type="spearman")
# pull out correlation values
male.socF.R <- male.socF.table$r
# pull out p-vals
male.socF.p <- male.socF.table$P

## Define notions for significance levels; spacing is important.
mystars <- ifelse(male.socF.p < .0001, "****", 
                  ifelse(male.socF.p < .001, "*** ",
                         ifelse(male.socF.p < .01, "**  ",
                                ifelse(male.socF.p < .05, "*    ", "    "))))

## truncate the correlation matrix to two decimal
male.socF.R2 <- format(round(cbind(rep(-1.11, ncol(male.socF.cor)), 
                                   male.socF.R), 2))[,-1]

## build a new matrix that includes the correlations with their appropriate stars
Rnew <- matrix(paste(male.socF.R2, mystars, sep=""), ncol=ncol(male.socF.cor))
diag(Rnew) <- paste(diag(male.socF.R), " ", sep="")
rownames(Rnew) <- colnames(male.socF.cor)
colnames(Rnew) <- paste(colnames(male.socF.cor), "", sep="")

# Hide upper triangle
male.socF.upper<-Rnew
male.socF.upper[upper.tri(Rnew)]<-""
male.socF.upper<-as.data.frame(male.socF.upper)

write.csv(male.socF.upper, "output-files/male socF corr table_2025-03-13.csv")

#-------------------------------------------------------------------------------
# Summary stats tables
#-------------------------------------------------------------------------------

## all females

# filter for numeric columns
fem.numeric <- fem[ ,c(4,31,8:11,36,13,18,20,38,39,25,30,21:24)]
# make summary stat table (empty)
summary.stats.fem <- as.data.frame(matrix(NA, nrow=3, ncol=18))
colnames(summary.stats.fem) <- colnames(fem.numeric)
# sample size
summary.stats.fem[1,] <- colSums(!is.na(fem.numeric))
# mean
summary.stats.fem[2,] <- apply(fem.numeric, 2 ,mean,na.rm=T)
# sd
summary.stats.fem[3,] <- apply(fem.numeric, 2, sd, na.rm=T)
# transpose
t.summary.stats.fem <- as.data.frame(t(as.matrix(summary.stats.fem)))
# add column labels
colnames(t.summary.stats.fem) <- c("N","mean","sd")
# export table
write.csv(t.summary.stats.fem, "output-files/summary stats for all females_2025-03-13.csv")


### All males

# filter for numeric columns
male.numeric <- male[ ,c(4,32,8:11,14,16,37,39,21,31,17:20,36)]
# make summary stat table (empty)
summary.stats.male <- as.data.frame(matrix(NA, nrow=3, ncol=17))
colnames(summary.stats.male) <- colnames(male.numeric)
# sample size
summary.stats.male[1,] <- colSums(!is.na(male.numeric))
# mean
summary.stats.male[2,] <- apply(male.numeric,2,mean,na.rm=T)
# sd
summary.stats.male[3,] <- apply(male.numeric, 2, sd, na.rm=T)
# transpose
t.summary.stats.male <- as.data.frame(t(as.matrix(summary.stats.male)))
# add column labels
colnames(t.summary.stats.male) <- c("N","mean","sd")
# export table
write.csv(t.summary.stats.male, "output-files/summary stats for all males_2025-03-13.csv")







