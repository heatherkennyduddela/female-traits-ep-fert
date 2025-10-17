
################################################################################
# Final models for female social nest success (fecundity)
# Heather Kenny-Duddela
# Oct 21, 2024
################################################################################

# libraries
library(dplyr)
library(ggplot2)
library(lubridate)
library(AICcmodavg) # for calculating AICc
library(car) # for Anova and vif
library(ggpubr) # for making plots with multiple panels

# load data
fem.season <- read.csv("input-files/female_season_table_hatched.csv")


# add scaled variables
scaled <- as.data.frame(lapply(fem.season[,c(8:11,21:24,30:31,36)], scale))
colnames.scaled <- paste(colnames(scaled),"scaled",sep="_")
colnames(scaled) <- colnames.scaled

fem <- cbind(fem.season, scaled)


#-------------------------------------------------------------------------------
# Check null models
#-------------------------------------------------------------------------------

# histogram of female fecundity
ggplot(fem, aes(x=tot_eggs)) + geom_histogram()

# check mean and variance
mean(fem$tot_eggs) #8.255319
var(fem$tot_eggs) #5.06383 reasonably close, though under dispersed

# compare null Poisson to null linear
# linear
null.lm <- lm(tot_eggs ~ 1, data=fem)
summary(null.lm)

# poisson
null.pois <- glm(tot_eggs ~ 1, data=fem, family="poisson")
summary(null.pois)

# check AICc
AICc(null.lm) #212.8819 simple linear is better!
AICc(null.pois) #216.9712

# check residuals
hist(null.lm$residuals)
plot(null.lm$fitted.values, null.lm$residuals)



#-------------------------------------------------------------------------------
# Fit models for female fecundity
#-------------------------------------------------------------------------------


## Fit the following models

# 1) female tail, male tail, female throat, male throat, female breast, male breast

### Mod 1
eggs.adj1 <- lm(tot_eggs ~ tail_scaled +
                  throat_avg_bright_scaled +
                  breast_avg_bright_scaled +
                  ci_1_julian_scaled, data=fem)

# get table of coefficients
table.eggs.adj1 <- summary(eggs.adj1)
table2.eggs.adj1 <- as.data.frame(cbind(table.eggs.adj1$coefficients, 
                              confint(eggs.adj1)))

# format CI as text
table2.eggs.adj1$ci.text <- paste("(", round(table2.eggs.adj1$`2.5 %`, 2), ", ", 
                                   round(table2.eggs.adj1$`97.5 %`, 2),")", sep="")

write.csv(table2.eggs.adj1, "output-files/fecundity final mod table scaled.csv")



### Mod 3, add fixed effect for site size
eggs.adj3 <- lm(tot_eggs ~ site.size + tail + socM_tail +
                  throat_avg_bright + socM_t.avg.bright +
                  breast_avg_bright + socM_r.avg.bright +
                  ci_1_julian, data=fem)

# get table of coefficients
table.eggs.adj3 <- summary(eggs.adj3)
# get table of vif values
vif.eggs.adj3 <- as.data.frame(vif(eggs.adj3))
vif2.eggs.adj3 <- rbind(NA, NA, vif.eggs.adj3)
# combine coefficients, CIs, and vif
table2.eggs.adj3 <- cbind(table.eggs.adj3$coefficients, 
                          confint(eggs.adj3), 
                          vif2.eggs.adj3)

write.csv(table2.eggs.adj3, "output-files/fecundity final mod table site.csv")


## Mod 4, no CI
eggs.adj4 <- lm(tot_eggs ~ tail + socM_tail +
                  throat_avg_bright + socM_t.avg.bright +
                  breast_avg_bright + socM_r.avg.bright,
                data=fem)

# get table of coefficients
table.eggs.adj4 <- summary(eggs.adj4)
# get table of vif values
vif.eggs.adj4 <- as.data.frame(vif(eggs.adj4))
vif2.eggs.adj4 <- rbind(NA, vif.eggs.adj4)
# combine coefficients, CIs, and vif
table2.eggs.adj4 <- cbind(table.eggs.adj4$coefficients, 
                          confint(eggs.adj4), 
                          vif2.eggs.adj4)

write.csv(table2.eggs.adj4, "output-files/fecundity final mod table no ci.csv")

## Mod 5, no male traits or CI
eggs.adj5 <- lm(tot_eggs ~ tail +
                  throat_avg_bright  +
                  breast_avg_bright,
                data=fem)

# get table of coefficients
table.eggs.adj5 <- summary(eggs.adj5)
# get table of vif values
vif.eggs.adj5 <- as.data.frame(vif(eggs.adj5))
vif2.eggs.adj5 <- rbind(NA, vif.eggs.adj5)
# combine coefficients, CIs, and vif
table2.eggs.adj5 <- cbind(table.eggs.adj5$coefficients, 
                          confint(eggs.adj5), 
                          vif2.eggs.adj5)

write.csv(table2.eggs.adj5, "output-files/fecundity final mod table no male traits.csv")



#-------------------------------------------------------------------------------
### Visualize results
## NOTE: run function for calculating CIs for glm


## female tail------------------------------------------------------------------

# range of tail values
pred.eggs.tail <- custom.newdata(old.data=fem, data.type="female", scaled=T,
                                 xvar="tail_scaled", xunscale="tail")

# add CIs and fit
eggs.tail.fit.ci <- predict(eggs.adj1, newdata = pred.eggs.tail[,1:8],
                type="response", interval="confidence")

# combine
pred.eggs.tail.adj1 <- cbind(pred.eggs.tail, eggs.tail.fit.ci)

# plot
eggs.tail <- ggplot(pred.eggs.tail.adj1) +
  geom_line(aes(x=tail, y=fit), size=1.5, color="#F2AA84", linetype="dashed") +
  geom_ribbon(aes(x=tail, ymin=lwr, ymax=upr), fill="#F2AA84",
              alpha=0.3) +
  geom_point(data=fem, aes(x=tail, y=tot_eggs), 
             alpha=0.5, color="#F2AA84") +
  xlab("Tail streamer length (mm)") +
  ylab("Total eggs laid") +
  theme_light() +
  ggtitle("    ")

ggsave("output-files/fecundity pred fem tail.png", h=3.5, w=4, scale=0.5)


## female breast----------------------------------------------------------------

# new data
pred.eggs.breast <- custom.newdata(old.data=fem, data.type="female", scaled=T,
                                 xvar="breast_avg_bright_scaled", 
                                 xunscale="breast_avg_bright")

# add CIs and fit
eggs.breast.fit.ci <- predict(eggs.adj1, newdata = pred.eggs.breast[,1:8],
                            type="response", interval="confidence")

# combine
pred.eggs.breast.adj1 <- cbind(pred.eggs.breast, eggs.breast.fit.ci)

# plot
eggs.breast <- ggplot(pred.eggs.breast.adj1) +
  geom_line(aes(x=breast_avg_bright, y=fit), size=1.5, color="#F2AA84",
            linetype="dashed") +
  geom_ribbon(aes(x=breast_avg_bright, ymin=lwr, ymax=upr), fill="#F2AA84",
              alpha=0.3) +
  geom_point(data=fem, aes(x=breast_avg_bright, y=tot_eggs), 
             alpha=0.5, color="#F2AA84") +
  xlab("Breast average brightness") +
  ylab(NULL) +
  theme_light() +
  ggtitle("    ")

ggsave("output-files/fecundity pred fem breast.png", h=3.5, w=4, scale=0.5)


## female throat ---------------------------------------------------------------

# new data
# range of tail values
pred.eggs.throat <- custom.newdata(old.data=fem, data.type="female", scaled=T,
                                 xvar="throat_avg_bright_scaled", 
                                 xunscale="throat_avg_bright")


# add CIs and fit
eggs.throat.fit.ci <- predict(eggs.adj1, newdata = pred.eggs.throat[,1:8],
                            type="response", interval="confidence")

# combine
pred.eggs.throat.adj1 <- cbind(pred.eggs.throat, eggs.throat.fit.ci)

# plot
eggs.throat <- ggplot(pred.eggs.throat.adj1) +
  geom_line(aes(x=throat_avg_bright, y=fit), size=1.5, 
            color="#F2AA84", linetype="dashed") +
  geom_ribbon(aes(x=throat_avg_bright, ymin=lwr, ymax=upr), 
              fill="#F2AA84",alpha=0.3) +
  geom_point(data=fem, aes(x=throat_avg_bright, y=tot_eggs), 
             alpha=0.5, color="#F2AA84") +
  xlab("Throat average brightness") +
  ylab(NULL) +
  theme_light() +
  ggtitle("    ")

ggsave("output-files/fecundity pred fem throat.png", h=3.5, w=4, scale=0.5)

## plot the three graphs together
ggarrange(eggs.tail, eggs.throat, eggs.breast, nrow=1, ncol=3, align="hv",
          labels=c("a","b","c"),label.x=0.2, label.y=1)

ggsave("output-files/effects on fecundity plots_2025-03-13_tan.png", w=7, h=2.5)

################################################################################
# correlations between fecundity and EPP
################################################################################

## binary EPP and fecundity
ggplot(fem, aes(y=tot_eggs, x=as.factor(ep.yes))) +
  geom_boxplot() + geom_point(position=position_jitter(h=0, w=0.1), 
                              alpha=0.3)

t.test(fem$tot_eggs ~ fem$ep.yes)

## Number EPO and fecundity
ggplot(fem, aes(x=tot_eggs, y=num.epo)) +
  geom_point(alpha=0.3, position=position_jitter(h=0.1, w=0.1)) +
  geom_smooth(method=lm, se=F)

cor.test(fem$tot_eggs, fem$num.epo, method="spearman")

## Proportion EPO and fecundity
ggplot(fem, aes(x=tot_eggs, y=season.prop.ep)) +
  geom_point(alpha=0.3) +
  geom_smooth(method=lm, se=F)

cor.test(fem$tot_eggs, fem$season.prop.ep, method="spearman")

## Number of mates
ggplot(fem, aes(x=tot_eggs, y=num.tot.mates)) +
  geom_point(position=position_jitter(h=0.1, w=0), alpha=0.3) +
  geom_smooth(method=lm, se=F)

cor.test(fem$tot_eggs, fem$num.tot.mates, method="spearman")



