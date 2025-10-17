
################################################################################
# Exploratory figures for female seasonal reproductive success
# Heather Kenny-Duddela
# Aug 20, 2024
################################################################################

# load data
fem <- read.csv("input-files/female_season_table_2025-03-13.csv")

male <- read.csv("input-files/male_season_table_2025-03-13.csv")

### libraries
library(tidyverse) # ggplot2, dplyr, lubridate
library(ggpubr) # for arranging multiple plots in the same panel

# add julian date for ci
fem$ci_1 <- ymd(fem$ci_1)
fem$ci_1_julian <- yday(fem$ci_1)

male$socF_ci <- ymd(male$socF_ci)
male$socF_ci_julian <- yday(male$socF_ci)



################################################################################

fem$ep.yes.factor <- factor(fem$ep.yes, levels=c(0,1), labels=c("No","Yes"))

### Do females with EPP fledge more chicks? 
fem.epp.rs <- ggplot(fem, aes(x=ep.yes.factor, y=tot_fledge, group=ep.yes)) +
  geom_boxplot(fill="#F2AA84") +
  geom_point(alpha=0.5,
             position=position_jitter(width=0.2, height=0)) +
  xlab("Female EP status") +
  ylab("Total offspring fledged") +
  ggtitle("t-test: t=-2.21, df=18.96, p=0.039") +
  ylim(0,10) +
  theme_light()

# ggsave("output-files/female RS and binary EPP.png", w=4, h=6, scale=0.5)
  

t.test(tot_fledge ~ ep.yes, data=fem)
# Welch Two Sample t-test
# 
# data:  tot_fledge by ep.yes
# t = -2.2141, df = 18.958, p-value = 0.03928
# alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
# 95 percent confidence interval:
#   -3.30994821 -0.09276673
# sample estimates:
#   mean in group 0 mean in group 1 
# 4.769231        6.470588 

cor.test(fem$tot_fledge, fem$ep.yes, method="spearman")
# Spearman's rank correlation rho
# 
# data:  fem$tot_fledge and fem$ep.yes
# S = 12207, p-value = 0.04471
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2942075 


### Do females with higher proportion EPO fledge more chicks? 
fem.prop.rs <- ggplot(fem, aes(x=season.prop.ep, y=tot_fledge, group=ep.yes)) +
  geom_point(alpha=0.5, 
             position=position_jitter(width=0.0, height=0.1)) +
  xlab("Female proportion EPO") +
  ylab("Total offspring fledged")+
  ggtitle("Spearman Rho=0.03, p=0.842") +
  geom_smooth(method=lm, se=F, linetype="dashed", color="#F2AA84")+
  theme_light() 

# ggsave("output-files/female rs and prop EPO.png", h=4, w=6.6, scale=0.5)

cor.test(fem$season.prop.ep, fem$tot_fledge, method="spearman")
# Spearman's rank correlation rho
# 
# data:  fem$season.prop.ep and fem$tot_fledge
# S = 16781, p-value = 0.8424
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.02980223 

## check spline smooth
ggplot(fem, aes(x=season.prop.ep, y=tot_fledge, group=ep.yes)) +
  geom_point(alpha=0.5, 
             position=position_jitter(width=0.0, height=0.1)) +
  xlab("Proportion EPO across the season") +
  ylab("Number of offspring fledged") +
  ggtitle("Female total RS and proportion EPO\n(R=0.16, p=0.245)") +
  geom_smooth(method=loess)



### Do females with a higher number of EPO fledge more chicks? 
fem.num.rs <- ggplot(fem, aes(x=num.epo, y=tot_fledge)) +
  geom_point(alpha=0.5, position=position_jitter(h=0.1, w=0.1)) +
  geom_smooth(method=lm, se=F, color="#F2AA84") +
  ylab(NULL) +
  xlab("Number of EP offspring") +
  ylab("Total offspring fledged") +
  ggtitle("Spearman Rho=0.30, p=0.039") +
  theme_light()

# ggsave("output-files/female RS and number EPO.png", h=4, w=6.6, scale=0.5)

cor.test(fem$num.epo, fem$tot_fledge, method="spearman")
# Spearman's rank correlation rho
# 
# data:  fem$num.epo and fem$tot_fledge
# S = 12065, p-value = 0.03879
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#      rho 
# 0.302467 


### Do female who lay earlier fledge more chicks?

ggplot(fem, aes(x=ci_1, y=tot_fledge)) +
         geom_point(alpha=0.5) +
  geom_smooth(method=lm) +
  xlab("First clutch lay date") +
  ylab("Total chicks fledged") +
  ggtitle("Female total RS vs. lay date\n(R=-0.53, p<0.001)")

ggsave("output-files/female RS vs lay date.png")

cor.test(fem$ci_1_julian, fem$tot_fledge, method="spearman")
# Spearman's rank correlation rho
# 
# data:  fem$ci_1_julian and fem$tot_fledge
# S = 26436, p-value = 0.0001345
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.5284522 

# check relationship with July dates removed
fem.early <- subset(fem, fem$ci_1_julian<187)

ggplot(fem.early, aes(x=ci_1, y=tot_fledge)) +
  geom_point(alpha=0.5) +
  geom_smooth(method=lm)


### do older females have more offspring? 

# subset known age females
fem.age <- subset(fem, fem$certainty!="low")

ggplot(fem.age, aes(x=age.class, y=tot.chick)) +
  geom_boxplot() +
  geom_point(alpha=0.5, 
             position=position_jitter(width=0.2, height=0)) +
  xlab("Female age class") +
  ylab("Total chicks fledged") +
  ggtitle("Female total RS vs age class\n(t-test p<0.001)")

ggsave("output-files/female RS vs age class.png")

t.test(tot.chick ~ age.class, data=fem.age)  
# Welch Two Sample t-test
# 
# data:  tot.chick by age.class
# t = 4.0067, df = 19.383, p-value = 0.0007296
# alternative hypothesis: true difference in means between group older and group young is not equal to 0
# 95 percent confidence interval:
#   1.614017 5.134646
# sample estimates:
#   mean in group older mean in group young 
# 7.727273            4.352941 


## Do females who lay more eggs fledge more offspring? 

fem.eggs.rs <- ggplot(fem, aes(x=tot_eggs, y=tot_fledge)) +
  geom_point(alpha=0.5, position=position_jitter(h=0.05, w=0.05)) +
  geom_smooth(method=lm, se=F, color="#F2AA84")+
  ylab("Total offspring fledged") +
  xlab("Total eggs laid") +
  ggtitle("Spearman Rho=0.47, p=0.004") +
  theme_light()

# ggsave("output-files/female RS and fecundity.png", h=6, w=8, scale=0.5)

cor.test(fem$tot_eggs, fem$tot_fledge, method="spearman")
# Spearman's rank correlation rho
# 
# data:  fem$tot_eggs and fem$tot_fledge
# S = 10225, p-value = 0.004326
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.4088508

ggarrange(fem.epp.rs, fem.num.rs, fem.prop.rs, fem.eggs.rs, 
          nrow=2, ncol=2)

## Do females who have more mates fledge more offspring?
fem.mates.rs <- ggplot(fem, aes(x=num.tot.mates, y=tot_fledge)) +
  geom_point(alpha=0.5, position=position_jitter(h=0.05, w=0.05))+
  geom_smooth(method=lm, se=F, color="#F2AA84")+
  ylab("Total offspring fledged") +
  xlab("Total genetic sires") +
  ggtitle("Spearman Rho=0.39, p=0.007") +
  theme_light()

cor.test(fem$num.tot.mates, fem$tot_fledge, method="spearman")
# Spearman's rank correlation rho
# 
# data:  fem$num.tot.mates and fem$tot_fledge
# S = 10636, p-value = 0.007523
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3850706 


##################################################################################
# Exploratory figs for male RS

# Note that for males, we estimate RS based on the sampled chicks only, not 
# on the number that actually fledged. 

male$ep.yes.factor <- factor(male$ep.yes, levels=c(0,1), labels=c("No","Yes"))

### Do males with EP sire more chicks overall? 
male.ep.rs <- ggplot(male, aes(x=ep.yes.factor, y=tot.chick)) +
  geom_boxplot(fill="#0033CC") + ylim(0,16) +
  geom_point(alpha=0.5,
             position=position_jitter(width=0.2, height=0))+
  xlab("Male EP status (outside social nest)") +
  ylab("Total offspring sired\n(EPO + WPO)") +
  ggtitle("t-test: t=-3.20, df=21.65, p=0.004")+
  theme_light()

# ggsave("output-files/male RS vs binary epp.png", w=4, h=6, scale=0.5)

t.test(tot.chick ~ ep.yes, data=male)
# Welch Two Sample t-test
# 
# data:  tot.chick by ep.yes
# t = -3.2047, df = 21.647, p-value = 0.004144
# alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
# 95 percent confidence interval:
#   -6.118842 -1.308103
# sample estimates:
#   mean in group 0 mean in group 1 
# 3.580645        7.294118 


### Do males with more EPO sire more chicks overall? 
male.num.rs <- ggplot(male, aes(x=ep.chick, y=tot.chick)) +
  geom_point(alpha=0.5, 
             position=position_jitter(h=0.1, w=0.1)) +
  geom_smooth(method=lm, se=F, color="#0033CC") +
  xlab("Number of EPO sired") +
  ylab("Total offspring sired\n(WPO + EPO)") +
  ggtitle("Spearman Rho=0.50, p<0.001") +
  theme_light()

# ggsave("output-files/male RS vs number EPO.png", h=4, w=6.6, scale=0.5)

cor.test(male$ep.chick, male$tot.chick, method="spearman")
# Spearman's rank correlation rho
# 
# data:  male$ep.chick and male$tot.chick
# S = 9187.5, p-value = 0.0002837
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5013314 


### Do males who sire a higher proportion of EP offspring have higher RS?
male.prop.rs <- ggplot(male, aes(x=prop.sired.epo, y=tot.chick)) +
  geom_point(alpha=0.5, position=position_jitter(h=0.1, w=0)) +
  geom_smooth(method=lm, se=F, linetype="solid", color="#0033CC") +
  xlab("Proportion of total sired that are EPO") +
  ylab("Total offspring sired\n(WPO + EPO)") +
  ggtitle("Spearman Rho=0.39, p=0.005") +
  theme_light()

# ggsave("output-files/male RS vs prop epo sired.png", h=4, w=6.6, scale=0.5)

cor.test(male$prop.sired.epo, male$tot.chick, method="spearman")
# Spearman's rank correlation rho
# 
# data:  male$prop.sired.epo and male$tot.chick
# S = 11169, p-value = 0.005624
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3937573 

## check spline smooth
ggplot(male, aes(x=prop.sired.epo, y=tot.chick)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="loess") +
  xlab("Proportion of EPO sired") +
  ylab("Total number of offspring sired (WPO + EPO)") +
  ggtitle("Male total RS vs proportion of EPO\n(R=0.18, p=0.22)")

# check pattern of tot chick vs. WPO
male$wp.chick <- male$tot.chick - male$ep.chick

ggplot(male, aes(x=wp.chick, y=tot.chick)) +
  geom_point(alpha=0.5) +
  geom_smooth(method=lm) +
  xlab("Number of WPO sired") +
  ylab("Total offspring sired (WPO + EPO)")+
  ggtitle("Male total RS vs number WPO\n(R=0.89, p<0.001)")

# ggsave("output-files/male RS vs number WPO.png")

cor.test(male$wp.chick, male$tot.chick, method="spearman")
# Spearman's rank correlation rho
# 
# data:  male$wp.chick and male$tot.chick
# S = 1236.7, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.9328739 


### Do males whose social mates lay earlier have higher RS?
ggplot(male, aes(x=ci_julian, y=tot.chick)) +
  geom_point(alpha=0.5) +
  geom_smooth(method=lm)+
  xlab("Social female first clutch lay date") +
  ylab("Total offspring sired (WPO + EPO)") +
  ggtitle("Male total RS vs social female lay date\n(R=-0.55, p<0.001)")

ggsave("output-files/male RS vs socF lay date.png")

cor.test(male$socF_ci_julian, male$tot.chick, method="spearman")
# Spearman's rank correlation rho
# 
# data:  male$socF_ci_julian and male$tot.chick
# S = 27485, p-value = 0.0003853
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.4917913 




### Do older males sire more offspring? 

# subset known age males
male.age <- subset(male, male$certainty!="low")

ggplot(male.age, aes(x=age.class, y=tot.chick)) +
  geom_boxplot() + ylim(0,16) +
  geom_point(alpha=0.5, 
             position=position_jitter(width=0.1, height=0)) +
  xlab("Male age class") +
  ylab("Total offspring sired (WPO + EPO)") +
  ggtitle("Male total RS vs age class for known age males\n(t-test p=0.009)")

ggsave("output-files/male RS vs age class.png")

t.test(tot.chick ~ age.class, data=male.age)
# Welch Two Sample t-test
# 
# data:  tot.chick by age.class
# t = 4.1437, df = 22.148, p-value = 0.0004199
# alternative hypothesis: true difference in means between group ASY and group SY is not equal to 0
# 95 percent confidence interval:
#   2.363293 7.095530
# sample estimates:
#   mean in group ASY  mean in group SY 
# 7.529412          2.800000 

##  Do males with more mates sire more offspring? 
male.mates.rs <- ggplot(male, aes(x=num.tot.mates, y=tot.chick)) +
  geom_point(alpha=0.5, position=position_jitter(h=0.1, w=0.1)) +
  geom_smooth(method=lm, se=F, linetype="solid", color="#0033CC") +
  xlab("Total genetic dams") +
  ylab("Total offspring sired\n(WPO + EPO)") +
  scale_x_continuous(breaks=c(1,2,3,4,5,6)) +
  ggtitle("Spearman Rho=0.63, p<0.001") +
  theme_light()

cor.test(male$num.tot.mates, male$tot.chick, method="spearman")
# Spearman's rank correlation rho
# 
# data:  male$num.tot.mates and male$tot.chick
# S = 6748.2, p-value = 1.332e-06
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.6337264 

# female and male RS in same panel ---------------------------------------------
#-------------------------------------------------------------------------------

ggarrange(fem.epp.rs, male.ep.rs, fem.num.rs, male.num.rs,
          fem.prop.rs, male.prop.rs, fem.mates.rs, male.mates.rs, fem.eggs.rs,
          nrow=5, ncol=2,
          labels=c("A","B","C","D","E","F","G","H","I"))

ggsave("output-files/RS patterns female and male combined_2025-10-07_bluetan.png", h=11, w=8)


