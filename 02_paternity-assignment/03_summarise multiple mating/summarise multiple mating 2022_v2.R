
#-------------------------------------------------------------------------------
# calculate number of mates for each female and male
#-------------------------------------------------------------------------------

# libraries
library(tidyverse)
library(ggplot2)

### load data tables

# fertilization types summarized by clutch ID
# from 02_summarise total RS and fertilization types
shared.fert.clutch2 <- read.csv("input-files/fert_2022_by_clutchID.csv")

# full kinship table with metadata
# from 01_paternity 2022 assignment with lcMLkin output
load("input-files/kin2_22.Rdata")

# list of kids with unknown dads
# from 01_paternity 2022 assignment with lcMLkin output
load("input-files/kids_unk_dads_2022.Rdata")

# families for paternity table, original input file
fam.22 <- read.csv("input-files/fam_clutch_2022_updated broods.csv")
# replace "-" in seqID with "_"
fam.22$Seq.ID.test <- gsub("-","_",fam.22$Seq.ID)



# count number of mates per male
mates.by.male <- shared.fert.clutch2 %>%
  group_by(Band_dad, FamilyID_dad, site_size_dad) %>%
  summarise(num.mates = length(unique(Band_mom)),
            num.clutches = length(unique(clutch_id_ind2)),
            num.chicks = sum(fert))

# subset to keep only EP kids, to count number of EP mates
shared.fert.clutch.ep <- subset(shared.fert.clutch2, 
                                shared.fert.clutch2$fert_type!="wp")

ep.mates.by.male <- shared.fert.clutch.ep %>%
  group_by(Band_dad, FamilyID_dad, site_size_dad) %>%
  summarise(num.ep.mates = length(unique(Band_mom)),
            num.ep.clutches = length(unique(clutch_id_ind2)),
            num.ep.chicks = sum(fert))

# calculate chicks per mate for males
mates.by.male$chicks.per.mate <- mates.by.male$num.chicks/mates.by.male$num.mates


#################################################################################
# for kids with unassigned dads, check whether they are full or half sibs

# take kids from the full kinship table that are in the unassigned dads table
dad.unk.sibs22 <- kin2.22[which(kin2.22$Ind2 %in% dad.unk22$Ind2),]

# filter to keep unassigned kids compared to unassigned kids
dad.unk.sibs22.2 <- dad.unk.sibs22[which(dad.unk.sibs22$Ind1 %in% dad.unk22$Ind2), ]

# filter to keep only pairs that are within the same family
# note that cases where only one kid within a family had an unknown dad are excluded
# from this table
dad.unk.sibs.fam22 <- subset(dad.unk.sibs22.2, dad.unk.sibs22.2$FamilyID_ind2 == 
                               dad.unk.sibs22.2$FamilyID_ind1)

## pull out cases where only one kid had an unknown dad within a family

# list of kids with unk dads that are in the sibs table
kid.sib.list <- c(dad.unk.sibs.fam22$Ind1, dad.unk.sibs.fam22$Ind2)

# all kids with unk dads (not just one with sibs)
kid.unk.list <- c(dad.unk.sibs22$Ind2)

# check which kids from the unk list are in the sibs list
kid.in.sib.list.index <- which(kid.unk.list %in% kid.sib.list)

# pull out kids from unk list that are NOT in the sibs list
kid.not.in.sib.list <- kid.unk.list[-kid.in.sib.list.index]

# find list of unique kids that are not in sibs list
# because it was only one kid in the clutch with an unknown dad
# 4 kids total
kid.not.in.sib.list.unique <- unique(kid.unk.list[-kid.in.sib.list.index])


### Loop to calculate number of unknown dads within each family (within each mom)


# object for number of families with unassigned kids
unk.fam22 <- as.data.frame(unique(dad.unk.sibs.fam22$FamilyID_ind1))
colnames(unk.fam22)[1] <- "FamilyID"
# add column to store number of dads
unk.fam22$num.dads <- NA

for (i in 1:length(unk.fam22$FamilyID)) {
  # subset by family id
  fam <- subset(dad.unk.sibs.fam22, dad.unk.sibs.fam22$FamilyID_ind1 == unk.fam22$FamilyID[i])
  
  # get relatedness data
  rel <- fam[,c(1,2,6)]
  
  # use reshape to go from a pairwise list to a distance matrix
  rel.r <- reshape(rel, direction="wide", idvar="Ind2", timevar="Ind1")
  
  # remove first column which is just labels
  rel.r.mat <- as.matrix(rel.r[,-1]) 
  
  # give full sibs as 0
  rel.r.bin <- rel.r.mat
  rel.r.bin[rel.r.mat >= 0.239] <- 0 # cutoff value from range of relatedness from 
  # genetic full sibs in 2022
  
  # give half sibs as 1
  rel.r.bin[rel.r.mat < 0.239] <- 1
  
  # calculate row products
  row.prod <- apply(rel.r.bin, 1, prod, na.rm=T)
  
  # sum up an add 1 to get the number of fathers
  num.dad <- sum(row.prod) + 1
  
  # save in storage
  unk.fam22$num.dads[i] <- num.dad
  
}

# add in kids not in the sib list

kid.not.in.sib.list.unique <- as.data.frame(kid.not.in.sib.list.unique)
colnames(kid.not.in.sib.list.unique)[1] <- "Seq.ID.test"

# add family ID
kid.not.in.sib.list.unique.fam <- left_join(kid.not.in.sib.list.unique, 
                                            fam.22[,c(5,11)], by="Seq.ID.test")
# add number of unk dads
kid.not.in.sib.list.unique.fam$num.dads <- 1

# combine two tables
unk.fam22.complete <- rbind(unk.fam22, kid.not.in.sib.list.unique.fam[,c(2:3)])


# count number of mates per female
mates.by.fem <- shared.fert.clutch2 %>%
  group_by(Band_mom, FamilyID_mom, site_size_mom) %>%
  summarise(num.mates = length(unique(Band_dad)),
            num.clutches = length(unique(clutch_id_ind2)),
            num.chicks = sum(fert))

### add in the unknown EP males

# loop to go through the UNK families and add to proper female
# the first UNK dad is already counted in the mates.by.fem table, so only 
# cases where there was more than one unknown need to be updated

mates.by.fem.complete <- mates.by.fem

for (i in 1:length(unk.fam22.complete$FamilyID)) {
  fam.index <- unk.fam22.complete$FamilyID[i]
  female.index <- which(mates.by.fem.complete$FamilyID_mom == fam.index)
  mates.sum <- (mates.by.fem.complete$num.mates[female.index] + unk.fam22.complete$num.dads[i] -1)
  mates.by.fem.complete$num.mates[female.index] <- mates.sum
}


### calculate number of EP mates per female

ep.mates.by.fem <- shared.fert.clutch.ep %>%
  group_by(Band_mom, FamilyID_mom, site_size_mom) %>%
  summarise(num.ep.mates = length(unique(Band_dad)),
            num.ep.clutches = length(unique(clutch_id_ind2)),
            num.ep.chicks = sum(fert))

# add in additional EP mates from unknown sires
# this time, just change the 2 families with multiple unknown
ep.mates.by.fem.complete <- ep.mates.by.fem

ep.mates.by.fem.complete$num.ep.mates[which(ep.mates.by.fem.complete$FamilyID_mom==
                                              "CHR-074")] <- 3
ep.mates.by.fem.complete$num.ep.mates[which(ep.mates.by.fem.complete$FamilyID_mom==
                                              "CHR-112")] <- 3

# calculate chicks per mate
mates.by.fem.complete$mates.per.chick <- mates.by.fem.complete$num.chicks/mates.by.fem.complete$num.mates

#-------------------------------------------------------------------------------
### Save tables with multiple mating

## Combine total mates and EP mates for males

mates.by.male.all <- left_join(mates.by.male, ep.mates.by.male, 
                               by=c("Band_dad","FamilyID_dad","site_size_dad"))
# for males that did not have any ep mates, fill in NAs with 0
mates.by.male.all$num.ep.mates[which(is.na(mates.by.male.all$num.ep.mates))] <- 0
mates.by.male.all$num.ep.clutches[which(
  is.na(mates.by.male.all$num.ep.clutches))] <- 0
mates.by.male.all$num.ep.chicks[which(is.na(mates.by.male.all$num.ep.chicks))] <- 0

# remove the NA row from unknown males
mates.by.male.all2 <- subset(mates.by.male.all, 
                             !is.na(mates.by.male.all$Band_dad))

# update column names
colnames(mates.by.male.all2)[4:7] <- c("num.tot.mates" ,"num.tot.clutches",
                                       "num.tot.chicks", "tot.chicks.per.mate")

# save table
write.csv(mates.by.male.all2, 
          "generated-files/male number of mates 2022.csv", row.names = F)

## combine total mates and EP mates for females

mates.by.fem.all <-left_join(mates.by.fem.complete, ep.mates.by.fem.complete,
                             by=c("Band_mom", "FamilyID_mom", "site_size_mom"))

# for females with no EP mates, fill in NAs with 0s
mates.by.fem.all$num.ep.mates[which(is.na(mates.by.fem.all$num.ep.mates))] <- 0
mates.by.fem.all$num.ep.clutches[which(
  is.na(mates.by.fem.all$num.ep.clutches))] <- 0
mates.by.fem.all$num.ep.chicks[which(is.na(mates.by.fem.all$num.ep.chicks))] <- 0

# remove the NA row from unsampled females
mates.by.fem.all2 <- subset(mates.by.fem.all, 
                             !is.na(mates.by.fem.all$Band_mom))

# update column names
colnames(mates.by.fem.all2)[4:6] <- c("num.tot.mates" ,"num.tot.clutches",
                                       "num.tot.chicks")

# save table
write.csv(mates.by.fem.all2, 
          "generated-files/female number of mates 2022.csv", row.names=F)


#-------------------------------------------------------------------------------
### Calculate number of mates per clutch (females only)

# initial count for all mates (ep and wp)
mates.by.fem.clutch <- shared.fert.clutch2 %>%
  group_by(Band_mom, FamilyID_mom, clutch_id_ind2, brood) %>%
  summarise(num.tot.mates = length(unique(Band_dad)),
            num.tot.chicks = sum(fert))

# just number of ep mates
# subset to keep only EP kids, to count number of EP mates
mates.by.fem.clutch.ep <- shared.fert.clutch.ep %>%
  group_by(Band_mom, FamilyID_mom, clutch_id_ind2, brood) %>%
  summarise(num.ep.mates = length(unique(Band_dad)),
            num.ep.chicks = sum(fert))


### Calculate number of unknown sires for each clutch

# filter to keep only pairs that are within the same clutch
# note that cases where only one kid within a clutch had an unknown dad are excluded
# from this table
dad.unk.sibs.clutch <- subset(dad.unk.sibs22.2, dad.unk.sibs22.2$clutch_id_ind2 == 
                               dad.unk.sibs22.2$clutch_id_ind1)

## pull out cases where only one kid had an unknown dad within a clutch

# list of kids with unk dads that are in the sibs table
kid.clutch.list <- c(dad.unk.sibs.clutch$Ind1, dad.unk.sibs.clutch$Ind2)

# all kids with unk dads (not just one with sibs)
kid.unk.list <- c(dad.unk.sibs22$Ind2)

# check which kids from the unk list are in the sibs list
kid.in.sib.list.clutch.index <- which(kid.unk.list %in% kid.clutch.list)

# pull out kids from unk list that are NOT in the sibs list
kid.not.in.sib.clutch.list <- kid.unk.list[-kid.in.sib.list.clutch.index]

# find list of unique kids that are not in sibs list
# because it was only one kid in the clutch with an unknown dad
# 8 kids total
kid.not.in.sib.clutch.list.unique <- unique(kid.unk.list[-kid.in.sib.list.clutch.index])

kid.not.in.sib.clutch.list.unique

### Loop to calculate number of unknown dads within each clutch


# object for number of clutches with unassigned kids
unk.clutch <- as.data.frame(unique(dad.unk.sibs.clutch$clutch_id_ind1))
colnames(unk.clutch)[1] <- "clutchID"
# add column to store number of dads
unk.clutch$num.dads <- NA

for (i in 1:length(unk.clutch$clutchID)) {
  # subset by clutch id
  fam <- subset(dad.unk.sibs.clutch, dad.unk.sibs.clutch$clutch_id_ind1 == unk.clutch$clutchID[i])
  
  # get relatedness data
  rel <- fam[,c(1,2,6)]
  
  # use reshape to go from a pairwise list to a distance matrix
  rel.r <- reshape(rel, direction="wide", idvar="Ind2", timevar="Ind1")
  
  # remove first column which is just labels
  rel.r.mat <- as.matrix(rel.r[,-1]) 
  
  # give full sibs as 0
  rel.r.bin <- rel.r.mat
  rel.r.bin[rel.r.mat >= 0.239] <- 0 # cutoff value from range of relatedness from 
  # genetic full sibs in 2022
  
  # give half sibs as 1
  rel.r.bin[rel.r.mat < 0.239] <- 1
  
  # calculate row products
  row.prod <- apply(rel.r.bin, 1, prod, na.rm=T)
  
  # sum up an add 1 to get the number of fathers
  num.dad <- sum(row.prod) + 1
  
  # save in storage
  unk.clutch$num.dads[i] <- num.dad
  
}

# There are no cases where more than one unknown sire was represented in each clutch.
# Therefore, we do not need to adjust the initial table of female mates per clutch

# combine total mates and just ep mates for the clutches
mates.by.fem.clutch.all <- left_join(mates.by.fem.clutch, mates.by.fem.clutch.ep, 
                                     by=c("Band_mom","FamilyID_mom",
                                          "clutch_id_ind2", "brood"))

# fill in NAs for females with zero EP mates
mates.by.fem.clutch.all$num.ep.mates[which(
  is.na(mates.by.fem.clutch.all$num.ep.mates))] <- 0

mates.by.fem.clutch.all$num.ep.chicks[which(
  is.na(mates.by.fem.clutch.all$num.ep.chicks))] <- 0

# remove the NA rows for the unsampled females
mates.by.fem.clutch.all2 <- mates.by.fem.clutch.all[-c(91,92), ]

# save table
write.csv(mates.by.fem.clutch.all2, row.names = F,
          "generated-files/female number of mates by clutch 2022.csv")


