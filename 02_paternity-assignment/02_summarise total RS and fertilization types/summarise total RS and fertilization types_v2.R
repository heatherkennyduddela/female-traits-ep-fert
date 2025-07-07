
################################################################################
# script for summarizing total RS and multiple mating by males and females
# Identifies genetic families and full sib relationships
# Makes plots of fertilization types based on sites size, brood, etc
# 
# Heather Kenny-Duddela
# July 14, 2023

# libraries
library(tidyverse)
library(ggplot2)

# load data - table of assigned parents
# from object in paternity 2022 assignment with lcMLkin output called kin.cutoff.info
load("input-files/kin_cutoff_info.Rdata")

# load full kinship table (not just parent-offspring relationships)
# from paternity 2022 assignment with lcMLkin output script
load("input-files/kin2_22.Rdata")


#-------------------------------------------------------------------------------
# Assess variation in relatedness value for full sibs
# Identify full sibs by having same mom and dad assigned
#-------------------------------------------------------------------------------


# switch to wide format where each row is a kid, and there are cols for dad and mom
# first just do dads
po.dad.wide <- subset(kin.cutoff.info, kin.cutoff.info$Type_ind1=="dad")
# pull out only cols that need to get "dad" labels added
po.dad.wide2 <- po.dad.wide[,c(1,3:11)]
colnames(po.dad.wide2)[2:10] <- c("k0_hat_dad", "k1_hat_dad", "k2_hat_dad", "pi_HAT_dad",
                                  "nbSNP_dad", "Band_dad","Site_dad", "Type_dad" , "FamilyID_dad")

# do same for moms
po.mom.wide <- subset(kin.cutoff.info, kin.cutoff.info$Type_ind1=="mom")
po.mom.wide2 <- po.mom.wide[,c(1,3:11)]
colnames(po.mom.wide2)[2:10] <- c("k0_hat_mom", "k1_hat_mom", "k2_hat_mom", "pi_HAT_mom",
                                  "nbSNP_mom", "Band_mom","Site_mom", "Type_mom" , "FamilyID_mom")

# pull out columns with relevant kid info, and no duplicated kids
po.kids <- kin.cutoff.info[-which(duplicated(kin.cutoff.info$Ind2)),c(1,15:21)]

# add dad columns to kid table
po.wide.d <- left_join(po.kids, po.dad.wide2, by="Ind2")

# add mom columns to kid table
po.wide <- left_join(po.wide.d, po.mom.wide2, by="Ind2")

# create column to specify genetic family, kids with same genetic mom and dad
po.wide$genetic_fam <- paste(po.wide$Band_dad, po.wide$Band_mom, sep="_")

# move genetic family label near the front of the table
po.wide <- po.wide[,c(1, 27, 2:26)]

# remove duplicate rows due to some males having multiple family IDs
# identify duplicated rows across kid ID, kid family, and dad band
po.wide.dup <- po.wide[duplicated(po.wide[,c(1,6,15)]),]
# get all rows with duplicated kids, not just ones identified by duplicated
po.wide.dup2 <- po.wide[(po.wide$Ind2 %in% po.wide.dup$Ind2), ]
# make subset of kid table with duplicated kids and dads removed
po.wide.sub <- po.wide[-(which(po.wide$Ind2 %in% po.wide.dup$Ind2)), ]
# Reduce the duplicated kids table so that each kid-dad pair has only a single row
# For most cases there is a FamilyID_dad that matches the kid FamilyID, so preferentially take
# these matching FamilyIDs instead of non-matching ones. These are true within-pair offspring. 
# For CO_83279 there is not a matching FamilyID_dad, so this is an EP kid from a dad that was
# attending two other nests. Here just take the first occurrence, which will not be marked as 
# T from duplicated()
po.wide.dup3 <- subset(po.wide.dup2, po.wide.dup2$Ind2=="CO_83279" &
                         !duplicated(po.wide.dup2[,c(1,6,15)]) # keep first row of CO_83279
                       |
                         po.wide.dup2$Ind2!="CO_83279") # Keep all other kids
# Keep other kid rows where kid FamilyID matches dad FamilyID, and also keep the one CO_83279 row
po.wide.dup4 <- po.wide.dup3[which(po.wide.dup3$FamilyID_ind2 == po.wide.dup3$FamilyID_dad |
                                     po.wide.dup3$Ind2=="CO_83279"), ]
# Remove duplicate rows where kid FamilyID matches dad FamilyID (because these rows are true duplicates)
po.wide.dup5 <- po.wide.dup4[-which(duplicated(po.wide.dup4[,c(1,6,15)])),]

# add keeper rows back to the main dataframe
po.wide.c <- rbind(po.wide.sub, po.wide.dup5)


################################################################################

# add column for categorical site size
po.wide.c$site_type <- NA

po.wide.c$site_type[which(po.wide.c$Site_ind2 == "Cathy's" |
                            po.wide.c$Site_ind2 == "Karen's" |
                            po.wide.c$Site_ind2 == "Dome House" |
                            po.wide.c$Site_ind2 == "Marte's" |
                            po.wide.c$Site_ind2 == "Mary Ann's" |
                            po.wide.c$Site_ind2 == "Speedwell")] <- "solitary"

po.wide.c$site_type[which(po.wide.c$Site_ind2 == "Urban Farm Girlz" |
                            po.wide.c$Site_ind2 == "McCauley" |
                            po.wide.c$Site_ind2 == "Struthers")] <- "small"

po.wide.c$site_type[which(po.wide.c$Site_ind2 == "Boyer" |
                            po.wide.c$Site_ind2 == "Blue Cloud" |
                            po.wide.c$Site_ind2 == "Make Believe" |
                            po.wide.c$Site_ind2 == "Cooks")] <- "medium"          

po.wide.c$site_type[which(po.wide.c$Site_ind2 == "CHR")] <- "large"

# reorder factors
po.wide.c$site_type <- factor(po.wide.c$site_type,
                              levels=c("solitary", "small","medium","large"))

# rename brood column
colnames(po.wide.c)[8] <- "brood"

# save wide file
write.csv(po.wide.c, file="generated-files/kin_2022_parent_offspring_assigned_wide.csv")


################################################################################


# count up number of offspring produced by each genetic family
shared.fert.list <- po.wide.c %>% 
  group_by(genetic_fam, Band_dad, Band_mom, 
           Site_dad, Site_mom, FamilyID_dad, FamilyID_mom, site_type) %>%
  summarise(fert = n())

hist(shared.fert.list$fert)

# add labels to specify fertilization type: within-pair (wp), extra-pair same site (ep_same), 
# extra-pair different site (ep_diff), unknown dad (dad_unk), unknown mom (mom_unk)
# based on matching mom and dad FamilyIDs and Site 

shared.fert.list$fert_type <- NA
shared.fert.list$fert_type[which(shared.fert.list$FamilyID_dad == 
                                   shared.fert.list$FamilyID_mom)] <- "wp"
shared.fert.list$fert_type[which(shared.fert.list$FamilyID_dad != shared.fert.list$FamilyID_mom &
                                   shared.fert.list$Site_dad == shared.fert.list$Site_mom)] <- "ep_same"
shared.fert.list$fert_type[which(shared.fert.list$FamilyID_dad != shared.fert.list$FamilyID_mom &
                                   shared.fert.list$Site_dad != shared.fert.list$Site_mom)] <- "ep_diff"
shared.fert.list$fert_type[which(is.na(shared.fert.list$FamilyID_mom))] <- "mom_unk"
shared.fert.list$fert_type[which(is.na(shared.fert.list$FamilyID_dad))] <- "dad_unk"


# save shared fertilizations list
write.csv(shared.fert.list, file="generated-files/fert_2022_by_genetic_family.csv", row.names=F)



#-------------------------------------------------------------------------------
# at the nest/clutch level, look at proportion of outside site dads
# get a sense for whether patterns are usually lots of EP at few nests, or 
# some EP at most nests
#-------------------------------------------------------------------------------


# summarise at clutch level
shared.fert.clutch <- po.wide.c %>% 
  group_by(clutch_id_ind2, Band_dad, Band_mom, 
           Site_dad, Site_mom, FamilyID_dad, FamilyID_mom, brood, site_type) %>%
  summarise(fert = n())

# number of clutches
length(unique(shared.fert.clutch$clutch_id_ind2))

# assign fert categories
shared.fert.clutch$fert_type <- NA
shared.fert.clutch$fert_type[which(shared.fert.clutch$FamilyID_dad == 
                                     shared.fert.clutch$FamilyID_mom)] <- "wp"
shared.fert.clutch$fert_type[which(shared.fert.clutch$FamilyID_dad != shared.fert.clutch$FamilyID_mom &
                                     shared.fert.clutch$Site_dad == shared.fert.clutch$Site_mom)] <- "ep_same"
shared.fert.clutch$fert_type[which(shared.fert.clutch$FamilyID_dad != shared.fert.clutch$FamilyID_mom &
                                     shared.fert.clutch$Site_dad != shared.fert.clutch$Site_mom)] <- "ep_diff"
shared.fert.clutch$fert_type[which(is.na(shared.fert.clutch$FamilyID_mom))] <- "mom_unk"
shared.fert.clutch$fert_type[which(is.na(shared.fert.clutch$FamilyID_dad))] <- "dad_unk"


# summarise number and proportion of each offspring type
total.fert.clutch <- sum(shared.fert.clutch$fert)

# number of females where we have both 1st and 2nd broods
fem.both.broods <- shared.fert.clutch %>%
  group_by(FamilyID_mom, brood) %>%
  summarise(chicks = n())

fem.both.brood2 <- subset(fem.both.broods, fem.both.broods$brood == "1" |
                            fem.both.broods$brood == "2")

fem.both.brood3 <- fem.both.brood2 %>%
  group_by(FamilyID_mom) %>%
  summarise(brood_count = n())

# 33 females had multiple broods
length(subset(fem.both.brood3$FamilyID_mom, fem.both.brood3$brood_count==2))  

# proportion of within-pair offspring per clutch
shared.fert.clutch2 <- shared.fert.clutch %>% 
  group_by(clutch_id_ind2) %>%
  mutate(clutch_size = sum(fert))

shared.fert.wp <- subset(shared.fert.clutch2, shared.fert.clutch2$fert_type=="wp")
shared.fert.wp <- shared.fert.wp %>%
  mutate(prop.wp = fert/clutch_size)

shared.fert.wp2 <- left_join(shared.fert.clutch2, shared.fert.wp[,c(1,13)], by="clutch_id_ind2")
shared.fert.wp2$prop.wp[which(is.na(shared.fert.wp2$prop.wp))] <- 0


#-------------------------------------------------------------------------------
# pull out full-sib relatedness values based on genetic fams

# get list of unique genetic families with mom and dad IDs known

known.gen.fam <- subset(po.wide, !is.na(po.wide$FamilyID_dad) & !is.na(po.wide$FamilyID_mom))
num.gen.fam <- unique(known.gen.fam$genetic_fam)

# loop through the different families

# storage for full genetic sibs
sib.storage <- kin2.22[1,]
sib.storage[,] <- NA
sib.storage$genetic_fam <- NA

for (i in 1:length(num.gen.fam)) {
  fam <- subset(known.gen.fam, known.gen.fam$genetic_fam == num.gen.fam[i])
  if (length(unique(fam$Ind2)) > 1) {
    gen.sib <- subset(kin2.22, kin2.22$Ind1 %in% fam$Ind2 &
                        kin2.22$Ind2 %in% fam$Ind2)
    gen.sib$genetic_fam <- num.gen.fam[i]
    sib.storage <- rbind(sib.storage, gen.sib)}
}

# remove row of NAs
sib.storage <- subset(sib.storage, !is.na(sib.storage$Ind1))

# plot relatedness values for full sibs
plot(sib.storage$k0_hat, sib.storage$pi_HAT)

ggplot(sib.storage, aes(k0_hat, pi_HAT)) + geom_point(shape=1) +
  ylab("kinship value (min=0.244, max=0.572)") +
  xlab("Est probability of unrelated (min=0.206, max=0.637)") +
  ggtitle("Range of kinship and k0-hat values for full sibs")

ggsave("generated-files/full sib kinship values.png")

# max and min cutoff values
range(sib.storage$k0_hat) # k0 from 0.206 to 0.637
range(sib.storage$pi_HAT) # pi_hat from 0.244 to 0.572

full.sib.pi.cutoff <- min(sib.storage$pi_HAT)

#-------------------------------------------------------------------------------
# add site sizes to shared.fert.clutch2
#-------------------------------------------------------------------------------

# Blue Cloud - 7
# Cooks - 9
# CHR - 32
# Make Believe - 7
# Struthers - 3
# McCauley - 2
# Urban Farm Girlz - 2
# Dome House - 1
# Marte's - 1
# Boyer - 5
# Cathy's - 1
# Mary Ann's - 1

# for dads
shared.fert.clutch2$site_size_dad <- NA
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Blue Cloud")] <- 7
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Starlight")] <- 3
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Jay's")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Boyer")] <- 5
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Cathy's")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="CHR")] <- 32
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Mayas")] <- 9
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Cooks")] <- 9
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Reinarz")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Dome House")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Karen's")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Make Believe")] <- 7
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Marte's")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Mary Ann's")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="McCauley")] <- 2
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Speedwell")] <- 1
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Struthers")] <- 3
shared.fert.clutch2$site_size_dad[which(shared.fert.clutch2$Site_dad=="Urban Farm Girlz")] <- 2

# for moms
shared.fert.clutch2$site_size_mom <- NA
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Blue Cloud")] <- 7
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Starlight")] <- 3
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Jay's")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Boyer")] <- 5
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Cathy's")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="CHR")] <- 32
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Mayas")] <- 9
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Cooks")] <- 9
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Reinarz")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Dome House")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Karen's")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Make Believe")] <- 7
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Marte's")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Mary Ann's")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="McCauley")] <- 2
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Speedwell")] <- 1
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Struthers")] <- 3
shared.fert.clutch2$site_size_mom[which(shared.fert.clutch2$Site_mom=="Urban Farm Girlz")] <- 2


####################################################################################
# Save key data tables

# save fertilization types by clutch ID
write.csv(shared.fert.clutch2, file="generated-files/fert_2022_by_clutchID.csv", row.names=F)








