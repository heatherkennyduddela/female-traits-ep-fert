
################################################################################
# Explore patterns between fecundity and EPP for females
# Heather Kenny-Duddela
# Feb 5, 2025
################################################################################

# load data
fem <- read.csv("output-files/female_season_table.csv")

# libraries
library(ggplot2)

#-------------------------------------------------------------------------------

# Scatterplots

## binary EPP and fecundity
ggplot(fem, aes(y=tot_eggs, x=as.factor(ep.yes))) +
  geom_boxplot() + geom_point(position=position_jitter(h=0, w=0.1), 
                              alpha=0.3)

t.test(fem$tot_eggs ~ fem$ep.yes)
# Welch Two Sample t-test
# 
# data:  fem$tot_eggs by fem$ep.yes
# t = 0.34924, df = 18.302, p-value = 0.7309
# alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
# 95 percent confidence interval:
#   -1.427798  1.997934
# sample estimates:
#   mean in group 0 mean in group 1 
# 8.461538        8.176471 


## Number EPO and fecundity
ggplot(fem, aes(x=tot_eggs, y=num.epo)) +
  geom_point(alpha=0.3, position=position_jitter(h=0.1, w=0.1)) +
  geom_smooth(method=lm, se=F)

cor.test(fem$tot_eggs, fem$num.epo, method="spearman")
# Spearman's rank correlation rho
# 
# data:  fem$tot_eggs and fem$num.epo
# S = 16122, p-value = 0.6504
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.0678599 


## Proportion EPO and fecundity
ggplot(fem, aes(x=tot_eggs, y=season.prop.ep)) +
  geom_point(alpha=0.3) +
  geom_smooth(method=lm, se=F)

cor.test(fem$tot_eggs, fem$season.prop.ep, method="spearman")
# Spearman's rank correlation rho
# 
# data:  fem$tot_eggs and fem$season.prop.ep
# S = 19386, p-value = 0.4184
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.1208507 


## Number of mates
ggplot(fem, aes(x=tot_eggs, y=num.tot.mates)) +
  geom_point(position=position_jitter(h=0.1, w=0), alpha=0.3) +
  geom_smooth(method=lm, se=F)

cor.test(fem$tot_eggs, fem$num.tot.mates, method="spearman")
# Spearman's rank correlation rho
# 
# data:  fem$tot_eggs and fem$num.tot.mates
# S = 15645, p-value = 0.5234
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.09542715 




