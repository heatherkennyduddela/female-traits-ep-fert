
################################################################################
# Script to explore patterns across sites
# Heather Kenny-Duddela
# Nov 12, 2024
################################################################################

# libraries
library(dplyr)
library(ggplot2)
library(stringr)

# load data
fem.s <- read.csv("output-files/female_season_table_2025-03-13.csv")
male.s <- read.csv("output-files/male_season_table_2025-03-13.csv")
group <- read.csv("input-files/Group sizes for 2022.csv")

# add numeric group sizes
# first standardize site names
group$site[2] <- "Colorado Horse Rescue"
fem <- left_join(fem.s, group[,1:2], by="site")

male <- left_join(male.s, group[,1:2], by="site")


# Histogram of females across group size
ggplot(fem, aes(x=brood1.size)) + geom_histogram(binwidth=1, fill="lightblue",
                                                 color="black") +
  ggtitle("Distribution of individual females across group sizes") +
  xlab("Number of breeding pairs")

# histogram of group sizes
group2 <- unique(select(fem,site,brood1.size))

ggplot(group2, aes(x=brood1.size)) + geom_histogram(binwidth=1, fill="lightblue",
                                                   color="black") +
  ggtitle("Distribution of group sizes across sites") +
  xlab("Number of breeding pairs")

ggsave("output-files/distribution of group sizes across sites.png", h=2, w=4)


# female tail ------------------------------------------------------------------

## by site
ggplot(fem, aes(x=tail, y=site)) + geom_boxplot() +
  geom_point(shape=1)

## by small/med/large
ggplot(fem, aes(x=tail)) + geom_histogram(fill="lightblue", color="black") +
  facet_grid(site.size~.)

ggplot(fem, aes(x=tail, y=site.size)) + geom_boxplot() +
  geom_point(shape=1)

f.tail.site <- kruskal.test(tail ~ site, fem)
# Kruskal-Wallis rank sum test
# 
# data:  tail by site
# Kruskal-Wallis chi-squared = 10.011, df = 9, p-value = 0.3496

f.tail.size <- kruskal.test(tail ~ site.size, fem)
# Kruskal-Wallis rank sum test
# 
# data:  tail by site.size
# Kruskal-Wallis chi-squared = 3.7318, df = 2, p-value = 0.1548


# female throat ----------------------------------------------------------------

## by site
ggplot(fem, aes(x=throat_avg_bright, y=site)) + geom_boxplot() +
  geom_point(shape=1)

f.throat.site <- kruskal.test(throat_avg_bright ~ site, fem)
# Kruskal-Wallis rank sum test
# 
# data:  throat_avg_bright by site
# Kruskal-Wallis chi-squared = 10.816, df = 9, p-value = 0.2886


## by size
ggplot(fem, aes(x=throat_avg_bright, y=site.size)) + geom_boxplot() +
  geom_point(shape=1)

f.throat.size <- kruskal.test(throat_avg_bright ~ site.size, fem)
# Kruskal-Wallis rank sum test
# 
# data:  throat_avg_bright by site.size
# Kruskal-Wallis chi-squared = 2.9265, df = 2, p-value = 0.2315


# female breast ----------------------------------------------------------------

## by site
ggplot(fem, aes(x=breast_avg_bright, y=site)) + geom_boxplot() +
  geom_point(shape=1)

f.breast.site <- kruskal.test(breast_avg_bright ~ site, fem)
# Kruskal-Wallis rank sum test
# 
# data:  breast_avg_bright by site
# Kruskal-Wallis chi-squared = 13.394, df = 9, p-value = 0.1456


## by size
ggplot(fem, aes(x=breast_avg_bright, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

f.breast.size <- kruskal.test(breast_avg_bright ~ site.size, fem)
# Kruskal-Wallis rank sum test
# 
# data:  breast_avg_bright by site.size
# Kruskal-Wallis chi-squared = 7.0421, df = 2, p-value = 0.02957


# female belly -----------------------------------------------------------------

## by site
ggplot(fem, aes(x=belly_avg_bright, y=site)) + geom_boxplot() +
  geom_point(shape=1)

f.belly.site <- kruskal.test(belly_avg_bright ~ site, fem)
# Kruskal-Wallis rank sum test
# 
# data:  belly_avg_bright by site
# Kruskal-Wallis chi-squared = 9.3141, df = 9, p-value = 0.4088

## by size
ggplot(fem, aes(x=belly_avg_bright, y=site.size)) + geom_boxplot() +
  geom_point(shape=1)

f.belly.size <- kruskal.test(belly_avg_bright ~ site.size, fem)
# Kruskal-Wallis rank sum test
# 
# data:  belly_avg_bright by site.size
# Kruskal-Wallis chi-squared = 4.4332, df = 2, p-value = 0.109


# female vent -----------------------------------------------------------------

## by site
ggplot(fem, aes(x=vent_avg_bright, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

f.vent.site <- kruskal.test(vent_avg_bright ~ site, fem)
# Kruskal-Wallis rank sum test
# 
# data:  vent_avg_bright by site
# Kruskal-Wallis chi-squared = 13.816, df = 9, p-value = 0.129

## by size
ggplot(fem, aes(x=vent_avg_bright, y=site.size)) + geom_boxplot() +
  geom_point(shape=1)

f.vent.size <- kruskal.test(vent_avg_bright ~ site.size, fem)
# Kruskal-Wallis rank sum test
# 
# data:  vent_avg_bright by site.size
# Kruskal-Wallis chi-squared = 5.2622, df = 2, p-value = 0.072


# social male tail -------------------------------------------------------------

## by site
ggplot(fem, aes(x=socM_tail, y=site)) + geom_boxplot() +
  geom_point(shape=1)

socM.tail.site <- kruskal.test(socM_tail ~ site, fem)
# Kruskal-Wallis rank sum test
# 
# data:  socM_tail by site
# Kruskal-Wallis chi-squared = 2.3831, df = 8, p-value = 0.967

## by size
ggplot(fem, aes(x=socM_tail, y=site.size)) + geom_boxplot() +
  geom_point(shape=1)

socM.tail.size <- kruskal.test(socM_tail ~ site.size, fem)
# Kruskal-Wallis rank sum test
# 
# data:  socM_tail by site.size
# Kruskal-Wallis chi-squared = 0.20454, df = 2, p-value = 0.9028


# social male throat -----------------------------------------------------------

## by site
ggplot(fem, aes(x=socM_t.avg.bright, y=site)) + geom_boxplot() +
  geom_point(shape=1)

socM.throat.site <- kruskal.test(socM_t.avg.bright ~ site, fem)
# Kruskal-Wallis rank sum test
# 
# data:  socM_t.avg.bright by site
# Kruskal-Wallis chi-squared = 9.387, df = 9, p-value = 0.4023

## by size
ggplot(fem, aes(x=socM_t.avg.bright, y=site.size)) + geom_boxplot() +
  geom_point(shape=1)

socM.throat.size <- kruskal.test(socM_t.avg.bright ~ site.size, fem)
# Kruskal-Wallis rank sum test
# 
# data:  socM_t.avg.bright by site.size
# Kruskal-Wallis chi-squared = 2.4895, df = 2, p-value = 0.288


# social male breast -----------------------------------------------------------

## by site
ggplot(fem, aes(x=socM_r.avg.bright, y=site)) + geom_boxplot() +
  geom_point(shape=1)

socM.breast.site <- kruskal.test(socM_r.avg.bright ~ site, fem)
# Kruskal-Wallis rank sum test
# 
# data:  socM_r.avg.bright by site
# Kruskal-Wallis chi-squared = 5.6987, df = 8, p-value = 0.6809

## by size
ggplot(fem, aes(x=socM_r.avg.bright, y=site.size)) + geom_boxplot() +
  geom_point(shape=1)

socM.breast.size <- kruskal.test(socM_r.avg.bright ~ site.size, fem)
# Kruskal-Wallis rank sum test
# 
# data:  socM_r.avg.bright by site.size
# Kruskal-Wallis chi-squared = 3.2568, df = 2, p-value = 0.1962


# social male belly ------------------------------------------------------------

## by site
ggplot(fem, aes(x=socM_b.avg.bright, y=site)) + geom_boxplot() +
  geom_point(shape=1)

socM.belly.site <- kruskal.test(socM_b.avg.bright ~ site, fem)
# Kruskal-Wallis rank sum test
# 
# data:  socM_b.avg.bright by site
# Kruskal-Wallis chi-squared = 5.4115, df = 9, p-value = 0.7971

## by size
ggplot(fem, aes(x=socM_b.avg.bright, y=site.size)) + geom_boxplot() +
  geom_point(shape=1)

socM.belly.size <- kruskal.test(socM_b.avg.bright ~ site.size, fem)
# Kruskal-Wallis rank sum test
# 
# data:  socM_b.avg.bright by site.size
# Kruskal-Wallis chi-squared = 0.30696, df = 2, p-value = 0.8577


# social male vent -------------------------------------------------------------

## by site
ggplot(fem, aes(x=socM_v.avg.bright, y=site)) + geom_boxplot() +
  geom_point(shape=1)

socM.vent.site <- kruskal.test(socM_v.avg.bright ~ site, fem)
# Kruskal-Wallis rank sum test
# 
# data:  socM_v.avg.bright by site
# Kruskal-Wallis chi-squared = 6.3133, df = 9, p-value = 0.7082

## by size
ggplot(fem, aes(x=socM_v.avg.bright, y=site.size)) + geom_boxplot() +
  geom_point(shape=1)

socM.vent.size <- kruskal.test(socM_v.avg.bright ~ site.size, fem)
# Kruskal-Wallis rank sum test
# 
# data:  socM_v.avg.bright by site.size
# Kruskal-Wallis chi-squared = 3.2729, df = 2, p-value = 0.1947


# lay date ---------------------------------------------------------------------

## by site
ggplot(fem, aes(x=ci_1_julian, y=site)) + geom_boxplot() +
  geom_point(shape=1)

f.ci.site <- kruskal.test(ci_1_julian ~ site, fem)
# Kruskal-Wallis rank sum test
# 
# data:  ci_1_julian by site
# Kruskal-Wallis chi-squared = 9.7687, df = 9, p-value = 0.3695

## by size
ggplot(fem, aes(x=ci_1_julian, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

f.ci.size <- kruskal.test(ci_1_julian ~ site.size, fem)
# Kruskal-Wallis rank sum test
# 
# data:  ci_1_julian by site.size
# Kruskal-Wallis chi-squared = 1.4669, df = 2, p-value = 0.4803




# fecundity --------------------------------------------------------------------

## by site
ggplot(fem, aes(x=tot_eggs, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

eggs.site <- kruskal.test(tot_eggs ~ site, fem)
# Kruskal-Wallis rank sum test
# 
# data:  tot_eggs by site
# Kruskal-Wallis chi-squared = 9.0308, df = 9, p-value = 0.4344

## by size
ggplot(fem, aes(x=tot_eggs, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

eggs.site.size <- kruskal.test(tot_eggs ~ site.size, fem)
# Kruskal-Wallis rank sum test
# 
# data:  tot_eggs by site.size
# Kruskal-Wallis chi-squared = 2.1564, df = 2, p-value = 0.3402


# Proportion EPO ---------------------------------------------------------------

## by site
ggplot(fem, aes(x=season.prop.ep, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

f.propEP.site <- kruskal.test(season.prop.ep ~ site, fem)
# Kruskal-Wallis rank sum test
# 
# data:  season.prop.ep by site
# Kruskal-Wallis chi-squared = 7.4577, df = 9, p-value = 0.5896

## by size
ggplot(fem, aes(x=season.prop.ep, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

f.propEP.size <- kruskal.test(season.prop.ep ~ site.size, fem)
# Kruskal-Wallis rank sum test
# 
# data:  season.prop.ep by site.size
# Kruskal-Wallis chi-squared = 0.37749, df = 2, p-value = 0.828


# Number EPO -------------------------------------------------------------------

## by site
ggplot(fem, aes(x=num.epo, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

f.numEP.site <- kruskal.test(num.epo ~ site, fem)
# Kruskal-Wallis rank sum test
# 
# data:  num.epo by site
# Kruskal-Wallis chi-squared = 7.7201, df = 9, p-value = 0.5626

## by size
ggplot(fem, aes(x=num.epo, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

f.numEP.size <- kruskal.test(num.epo ~ site.size, fem)
# Kruskal-Wallis rank sum test
# 
# data:  num.epo by site.size
# Kruskal-Wallis chi-squared = 0.47058, df = 2, p-value = 0.7903


# Binary EPP -------------------------------------------------------------------

## by site
ggplot(fem, aes(y=site, fill=as.factor(ep.yes))) + 
  geom_bar(position="fill")

## by size
ggplot(fem, aes(y=site.size, fill=as.factor(ep.yes))) + 
  geom_bar(position="fill")


# Number of sires---------------------------------------------------------------

## by site
ggplot(fem, aes(x=num.tot.mates, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

f.mates.site <- kruskal.test(num.tot.mates ~ site, fem)


## by size
ggplot(fem, aes(x=num.tot.mates, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

f.mates.size <- kruskal.test(num.tot.mates ~ site.size, fem)



################################################################################
#-------------------------------------------------------------------------------
# Patterns for male table
#-------------------------------------------------------------------------------


# Male tail --------------------------------------------------------------------

## by site
ggplot(male, aes(x=tail, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

m.tail.site <- kruskal.test(tail ~ site, male)
# Kruskal-Wallis rank sum test
# 
# data:  tail by site
# Kruskal-Wallis chi-squared = 4.7084, df = 9, p-value = 0.859

## by size
ggplot(male, aes(x=tail, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

m.tail.size <- kruskal.test(tail ~ site.size, male)
# Kruskal-Wallis rank sum test
# 
# data:  tail by site.size
# Kruskal-Wallis chi-squared = 1.2382, df = 2, p-value = 0.5384


# Male throat ------------------------------------------------------------------

## by site
ggplot(male, aes(x=throat_avg_bright, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

m.throat.site <- kruskal.test(throat_avg_bright ~ site, male)
# Kruskal-Wallis rank sum test
# 
# data:  throat_avg_bright by site
# Kruskal-Wallis chi-squared = 9.7873, df = 10, p-value = 0.4593

## by size
ggplot(male, aes(x=throat_avg_bright, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

m.throat.size <- kruskal.test(throat_avg_bright ~ site.size, male)
# Kruskal-Wallis rank sum test
# 
# data:  throat_avg_bright by site.size
# Kruskal-Wallis chi-squared = 3.3939, df = 2, p-value = 0.1832


# Male breast ------------------------------------------------------------------

## by site
ggplot(male, aes(x=breast_avg_bright, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

m.breast.site <- kruskal.test(breast_avg_bright ~ site, male)
# Kruskal-Wallis rank sum test
# 
# data:  breast_avg_bright by site
# Kruskal-Wallis chi-squared = 6.4845, df = 9, p-value = 0.6906

## by size
ggplot(male, aes(x=breast_avg_bright, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

m.breast.size <- kruskal.test(breast_avg_bright ~ site.size, male)
# Kruskal-Wallis rank sum test
# 
# data:  breast_avg_bright by site.size
# Kruskal-Wallis chi-squared = 3.9728, df = 2, p-value = 0.1372


# Male belly -------------------------------------------------------------------

## by site
ggplot(male, aes(x=belly_avg_bright, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

m.belly.site <- kruskal.test(belly_avg_bright ~ site, male)
# Kruskal-Wallis rank sum test
# 
# data:  belly_avg_bright by site
# Kruskal-Wallis chi-squared = 6.3198, df = 10, p-value = 0.7877

## by size
ggplot(male, aes(x=belly_avg_bright, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

m.belly.size <- kruskal.test(belly_avg_bright ~ site.size, male)
# Kruskal-Wallis rank sum test
# 
# data:  belly_avg_bright by site.size
# Kruskal-Wallis chi-squared = 1.0553, df = 2, p-value = 0.59


# Male vent --------------------------------------------------------------------

## by site
ggplot(male, aes(x=vent_avg_bright, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

m.vent.site <- kruskal.test(vent_avg_bright ~ site, male)
# Kruskal-Wallis rank sum test
# 
# data:  vent_avg_bright by site
# Kruskal-Wallis chi-squared = 5.7859, df = 10, p-value = 0.8329

## by size
ggplot(male, aes(x=vent_avg_bright, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

m.vent.size <- kruskal.test(vent_avg_bright ~ site.size, male)
# Kruskal-Wallis rank sum test
# 
# data:  vent_avg_bright by site.size
# Kruskal-Wallis chi-squared = 2.7759, df = 2, p-value = 0.2496


# social female tail -----------------------------------------------------------

## by site
ggplot(male, aes(x=socF_tail, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

socF.tail.site <- kruskal.test(socF_tail ~ site, male)
# Kruskal-Wallis rank sum test
# 
# data:  socF_tail by site
# Kruskal-Wallis chi-squared = 10.611, df = 10, p-value = 0.3886

## by size
ggplot(male, aes(x=socF_tail, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

socF.tail.size <- kruskal.test(socF_tail ~ site.size, male)
# Kruskal-Wallis rank sum test
# 
# data:  socF_tail by site.size
# Kruskal-Wallis chi-squared = 3.9382, df = 2, p-value = 0.1396


# Social female throat ---------------------------------------------------------

## by site
ggplot(male, aes(x=socF_t.avg.bright, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

socF.throat.site <- kruskal.test(socF_t.avg.bright ~ site, male)
# Kruskal-Wallis rank sum test
# 
# data:  socF_t.avg.bright by site
# Kruskal-Wallis chi-squared = 13.799, df = 10, p-value = 0.1823

## by size
ggplot(male, aes(x=socF_t.avg.bright, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

socF.throat.size <- kruskal.test(socF_t.avg.bright ~ site.size, male)
# Kruskal-Wallis rank sum test
# 
# data:  socF_t.avg.bright by site.size
# Kruskal-Wallis chi-squared = 5.5949, df = 2, p-value = 0.06097


# Social female breast ---------------------------------------------------------

## by site
ggplot(male, aes(x=socF_r.avg.bright, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

socF.breast.site <- kruskal.test(socF_r.avg.bright ~ site, male)
# Kruskal-Wallis rank sum test
# 
# data:  socF_r.avg.bright by site
# Kruskal-Wallis chi-squared = 15.239, df = 10, p-value = 0.1236

## by size
ggplot(male, aes(x=socF_r.avg.bright, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

socF.breast.size <- kruskal.test(socF_r.avg.bright ~ site.size, male)
# Kruskal-Wallis rank sum test
# 
# data:  socF_r.avg.bright by site.size
# Kruskal-Wallis chi-squared = 7.1172, df = 2, p-value = 0.02848


# Social female belly ----------------------------------------------------------

## by site
ggplot(male, aes(x=socF_b.avg.bright, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

socF.belly.site <- kruskal.test(socF_b.avg.bright ~ site, male)
# Kruskal-Wallis rank sum test
# 
# data:  socF_b.avg.bright by site
# Kruskal-Wallis chi-squared = 11.403, df = 10, p-value = 0.327

## by size
ggplot(male, aes(x=socF_b.avg.bright, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

socF.belly.size <- kruskal.test(socF_b.avg.bright ~ site.size, male)
# Kruskal-Wallis rank sum test
# 
# data:  socF_b.avg.bright by site.size
# Kruskal-Wallis chi-squared = 2.9925, df = 2, p-value = 0.224


# social female vent -----------------------------------------------------------

## by site
ggplot(male, aes(x=socF_v.avg.bright, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

socF.vent.site <- kruskal.test(socF_v.avg.bright ~ site, male)
# Kruskal-Wallis rank sum test
# 
# data:  socF_v.avg.bright by site
# Kruskal-Wallis chi-squared = 14.476, df = 10, p-value = 0.1523

## by size
ggplot(male, aes(x=socF_v.avg.bright, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

socF.vent.size <- kruskal.test(socF_v.avg.bright ~ site.size, male)
# Kruskal-Wallis rank sum test
# 
# data:  socF_v.avg.bright by site.size
# Kruskal-Wallis chi-squared = 3.578, df = 2, p-value = 0.1671


# Proportion EPO ---------------------------------------------------------------

## by site
ggplot(male, aes(x=prop.sired.epo, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

m.propEP.site <- kruskal.test(prop.sired.epo ~ site, male)
# Kruskal-Wallis rank sum test
# 
# data:  prop.sired.epo by site
# Kruskal-Wallis chi-squared = 8.9433, df = 10, p-value = 0.5375

## by size
ggplot(male, aes(x=prop.sired.epo, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

m.propEP.size <- kruskal.test(prop.sired.epo ~ site.size, male)
# Kruskal-Wallis rank sum test
# 
# data:  prop.sired.epo by site.size
# Kruskal-Wallis chi-squared = 1.7907, df = 2, p-value = 0.4085


# Number EPO -------------------------------------------------------------------

## by site
ggplot(male, aes(x=ep.chick, y=site)) + geom_boxplot() +
  geom_point(alpha=0.5)

m.numEP.site <- kruskal.test(ep.chick ~ site, male)
# Kruskal-Wallis rank sum test
# 
# data:  ep.chick by site
# Kruskal-Wallis chi-squared = 6.9066, df = 10, p-value = 0.7342

## by size
ggplot(male, aes(x=ep.chick, y=site.size)) + geom_boxplot() +
  geom_point(alpha=0.5)

m.numEP.size <- kruskal.test(ep.chick ~ site.size, male)
# Kruskal-Wallis rank sum test
# 
# data:  ep.chick by site.size
# Kruskal-Wallis chi-squared = 1.9156, df = 2, p-value = 0.3837


# Binary EPP -------------------------------------------------------------------

## by site
ggplot(male, aes(y=site, fill=as.factor(ep.yes))) + 
  geom_bar(position="fill")

## by size
ggplot(male, aes(y=site.size, fill=as.factor(ep.yes))) + 
  geom_bar(position="fill")

# Number of dams ---------------------------------------------------------------

m.mates.site <- kruskal.test(num.tot.mates ~ site, male)

m.mates.size <- kruskal.test(num.tot.mates ~ site.size, male)


#-------------------------------------------------------------------------------
# Write results to a summary table

test.f <- list(f.tail.site, f.tail.size, f.throat.site, f.throat.size,  
                   f.breast.site, f.breast.size,f.belly.site, f.belly.size,
                   f.vent.site, f.vent.size,f.ci.site, f.ci.size, f.mates.site, 
                   f.mates.size, f.numEP.site, f.numEP.size, f.propEP.site, 
                   f.propEP.size)

test.socM <- list(socM.tail.site, socM.tail.size, socM.throat.site, socM.throat.size,  
                   socM.breast.site, socM.breast.size,socM.belly.site, socM.belly.size,
                   socM.vent.site, socM.vent.size)

test.m <- list(m.tail.site, m.tail.size, m.throat.site, m.throat.size,  
                   m.breast.site, m.breast.size,m.belly.site, m.belly.size,
                   m.vent.site, m.vent.size, m.mates.site, 
                   m.mates.size, m.numEP.site, m.numEP.size, m.propEP.site, 
                   m.propEP.size)

test.socF <- list(socF.tail.site, socF.tail.size, socF.throat.site, socF.throat.size,  
                   socF.breast.site, socF.breast.size,socF.belly.site, socF.belly.size,
                   socF.vent.site, socF.vent.size)

KW.table <- as.data.frame(matrix(nrow=54, ncol=4, NA))
colnames(KW.table) <- c("bird", "data.name", "Chi.sqr","p.val")

# table for female traits

for (i in 1:length(test.f)) {
  test <- test.f[[i]]
  KW.table$bird[i] <- "female"
  KW.table$data.name[i] <- test$data.name
  KW.table$Chi.sqr[i] <- test$statistic
  KW.table$p.val[i] <- test$p.value
}


# table for social males

for (i in 1:length(test.socM)) {
  test <- test.socM[[i]]
  KW.table$bird[i+18] <- "socM"
  KW.table$data.name[i+18] <- test$data.name
  KW.table$Chi.sqr[i+18] <- test$statistic
  KW.table$p.val[i+18] <- test$p.value
}


# table for male

for (i in 1:length(test.m)) {
  test <- test.m[[i]]
  KW.table$bird[i+28] <- "male"
  KW.table$data.name[i+28] <- test$data.name
  KW.table$Chi.sqr[i+28] <- test$statistic
  KW.table$p.val[i+28] <- test$p.value
}

# table for social females

for (i in 1:length(test.socF)) {
  test <- test.socF[[i]]
  KW.table$bird[i+44] <- "socF"
  KW.table$data.name[i+44] <- test$data.name
  KW.table$Chi.sqr[i+44] <- test$statistic
  KW.table$p.val[i+44] <- test$p.value
}

## split data.name into two columns
cols <- as.data.frame(str_split_fixed(KW.table$data.name, " by ", 2))
colnames(cols) <- c("variable", "grouping")

KW.table2 <- cbind(KW.table, cols)

write.csv(KW.table2, "output-files/Kruskal-Wallis results table.csv",
          row.names=F)

