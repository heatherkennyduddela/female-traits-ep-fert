################################################################################
# Script for fitting final models for proportion EPO
# Heather Kenny-Duddela
# Oct 16, 2024
################################################################################

### load data

# female table
fem.season <- read.csv("input-files/female_season_table_2025-03-13.csv")

# male table
male.season <- read.csv("input-files/male_season_table_2025-03-13.csv")

### libraries
library(tidyverse) #ggplot2 and dplyr
library(car) # for Anova and vif (variance inflation factor)
library(broom)# for calculating confidence intervals around coefficient estimates
library(AICcmodavg) # for calculating AIC
library(MASS) # for the glm.nb function

# NOTE: Also run the custom function scripts: 
# 1) function for custom new data.R
# 2) function to calculate glm CIs.R

#-------------------------------------------------------------------------------
# add scaled variables

scaled <- as.data.frame(lapply(fem.season[,c(8:11,21:24,30:31,36)], scale))
colnames.scaled <- paste(colnames(scaled),"scaled",sep="_")
colnames(scaled) <- colnames.scaled

fem <- cbind(fem.season, scaled)


scaled.m <- as.data.frame(lapply(male.season[,c(8:11,17:20,31:32,36)], scale))
colnames.scaled.m <- paste(colnames(scaled.m),"scaled",sep="_")
colnames(scaled.m) <- colnames.scaled.m

male <- cbind(male.season, scaled.m)


#-------------------------------------------------------------------------------
# Fit prop EPO models for females
#-------------------------------------------------------------------------------

# final model  

fem.prop.adj1 <- glm(season.prop.ep ~ tail_scaled + socM_tail_scaled +
                       throat_avg_bright_scaled + socM_t.avg.bright_scaled +
                       breast_avg_bright_scaled +socM_r.avg.bright_scaled +
                       ci_1_julian_scaled, 
                     data=fem, weights = tot.chick, family="binomial")

# get table of coefficients
table.fem.prop.adj1 <- summary(fem.prop.adj1)
table2.fem.prop.adj1 <- as.data.frame(cbind(table.fem.prop.adj1$coefficients, 
                                confint(fem.prop.adj1)))
# back transform
table2.fem.prop.adj1$beta.BT <- exp(table2.fem.prop.adj1$Estimate)
table2.fem.prop.adj1$lowerCI.BT <- exp(table2.fem.prop.adj1$`2.5 %`)
table2.fem.prop.adj1$upperCI.BT <- exp(table2.fem.prop.adj1$`97.5 %`)

write.csv(table2.fem.prop.adj1, "output-files/fem prop EPO final mod table scaled.csv")

# Add effect for site size
fem.prop.adj2 <- glm(season.prop.ep ~ site.size + tail_scaled + socM_tail_scaled +
                       throat_avg_bright_scaled + socM_t.avg.bright_scaled +
                       breast_avg_bright_scaled +socM_r.avg.bright_scaled +
                       ci_1_julian_scaled, 
                     data=fem, weights = tot.chick, family="binomial")

# get table of coefficients
table.fem.prop.adj2 <- summary(fem.prop.adj2)
table2.fem.prop.adj2 <- as.data.frame(cbind(table.fem.prop.adj2$coefficients, 
                              confint(fem.prop.adj2)))
# back transform
table2.fem.prop.adj2$beta.BT <- exp(table2.fem.prop.adj2$Estimate)
table2.fem.prop.adj2$lowerCI.BT <- exp(table2.fem.prop.adj2$`2.5 %`)
table2.fem.prop.adj2$upperCI.BT <- exp(table2.fem.prop.adj2$`97.5 %`)

write.csv(table2.fem.prop.adj2, 
          "output-files/fem prop EPO final mod table site scaled.csv")

## compare models with and without site
# support for keeping site in this model
anova(fem.prop.adj1, fem.prop.adj2)
# Analysis of Deviance Table
# 
# Model 1: season.prop.ep ~ tail_scaled + socM_tail_scaled + throat_avg_bright_scaled + 
#   socM_t.avg.bright_scaled + breast_avg_bright_scaled + socM_r.avg.bright_scaled + 
#   ci_1_julian_scaled
# Model 2: season.prop.ep ~ site.size + tail_scaled + socM_tail_scaled + 
#   throat_avg_bright_scaled + socM_t.avg.bright_scaled + breast_avg_bright_scaled + 
#   socM_r.avg.bright_scaled + ci_1_julian_scaled
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
# 1        35     99.413                       
# 2        33     93.182  2   6.2311  0.04435 *


# fit model without CI
fem.prop.adj3 <- glm(season.prop.ep ~ tail_scaled + socM_tail_scaled +
                       throat_avg_bright_scaled + socM_t.avg.bright_scaled +
                       breast_avg_bright_scaled +socM_r.avg.bright_scaled, 
                     data=fem, weights = tot.chick, family="binomial")
# get table of coefficients
table.fem.prop.adj3 <- summary(fem.prop.adj3)
table2.fem.prop.adj3 <- as.data.frame(cbind(table.fem.prop.adj3$coefficients, 
                              confint(fem.prop.adj3)))
# back transform
table2.fem.prop.adj3$beta.BT <- exp(table2.fem.prop.adj3$Estimate)
table2.fem.prop.adj3$lowerCI.BT <- exp(table2.fem.prop.adj3$`2.5 %`)
table2.fem.prop.adj3$upperCI.BT <- exp(table2.fem.prop.adj3$`97.5 %`)

write.csv(table2.fem.prop.adj3, "output-files/fem prop EPO table no ci scaled.csv")

###-----------------------------------------------------------------------------
# visualize results for female propEPO
###-----------------------------------------------------------------------------

## female tail------------------------------------------------------------------

# generate new data
pred.fem.prop.tail <- custom.newdata(old.data = fem, data.type="female", scaled=T,
                                     xvar="tail_scaled", 
                                     xunscale="tail")

# calculate fit and CIs
pred.fem.prop.tail.adj1 <- calculate.glm.ci(fem.prop.adj1, pred.fem.prop.tail)

ggplot(pred.fem.prop.tail.adj1) +
  geom_line(aes(x=tail, y=fit), size=1.5, linetype="solid",
            color="darkgreen") +
  geom_ribbon(aes(x=tail, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=tail, y=season.prop.ep), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Female tail length") +
  ylab("Seasonal proportion EPO")


# ggsave("output-files/fem prop EPO model pred fem tail.png", h=4, w=4.5)


## female throat----------------------------------------------------------------

# generate y-values
pred.fem.prop.throat <- custom.newdata(old.data = fem, data.type="female", scaled=T,
                                       xvar="throat_avg_bright_scaled", 
                                       xunscale="throat_avg_bright")

# calculate fit and CIs
pred.fem.prop.throat.adj1 <- calculate.glm.ci(fem.prop.adj1, 
                                              pred.fem.prop.throat)

ggplot(pred.fem.prop.throat.adj1) +
  geom_line(aes(x=throat_avg_bright, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=throat_avg_bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=throat_avg_bright, y=season.prop.ep), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Female throat average brightness") +
  ylab("Seasonal proportion EPO")


# ggsave("output-files/fem prop EPO model pred fem throat.png", h=4, w=4.5)


## female breast----------------------------------------------------------------

# generate y-values
pred.fem.prop.breast <- custom.newdata(old.data = fem, data.type="female", scaled=T,
                                       xvar="breast_avg_bright_scaled", 
                                       xunscale="breast_avg_bright")

# calcualte fit and CIs
pred.fem.prop.breast.adj1 <- calculate.glm.ci(fem.prop.adj1, 
                                              pred.fem.prop.breast)

# plot
ggplot(pred.fem.prop.breast.adj1) +
  geom_line(aes(x=breast_avg_bright, y=fit), size=1.5, linetype="solid",
            color="darkgreen") +
  geom_ribbon(aes(x=breast_avg_bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=breast_avg_bright, y=season.prop.ep), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Female breast average brightness") +
  ylab("Seasonal proportion EPO")

# ggsave("output-files/fem prop EPO model pred fem breast.png", h=4, w=4.5)


## social male tail-------------------------------------------------------------

# generate new data
pred.fem.prop.Mtail <- custom.newdata(old.data = fem, data.type="female", scaled=T,
                                      xvar="socM_tail_scaled", 
                                      xunscale="socM_tail")

# calculate fit and CIs
pred.fem.prop.Mtail.adj1 <- calculate.glm.ci(fem.prop.adj1, pred.fem.prop.Mtail)

ggplot(pred.fem.prop.Mtail.adj1) +
  geom_line(aes(x=socM_tail, y=fit), size=1.5, linetype="solid",
            color="darkgreen") +
  geom_ribbon(aes(x=socM_tail, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=socM_tail, y=season.prop.ep), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social male tail length") +
  ylab("Seasonal proportion EPO")

# ggsave("output-files/fem prop EPO model pred socM tail.png", h=4, w=4.5)


## social male throat-----------------------------------------------------------

# generate y-vals
pred.fem.prop.Mthroat <- custom.newdata(old.data = fem, data.type="female", scaled=T,
                                        xvar="socM_t.avg.bright_scaled", 
                                        xunscale="socM_t.avg.bright")

pred.fem.prop.Mthroat.adj1 <- calculate.glm.ci(fem.prop.adj1, pred.fem.prop.Mthroat)

ggplot(pred.fem.prop.Mthroat.adj1) +
  geom_line(aes(x=socM_t.avg.bright, y=fit), size=1.5,
            linetype="dashed") +
  geom_ribbon(aes(x=socM_t.avg.bright, ymax=upper, ymin=lower), alpha=0.3) +
  geom_point(data=fem, aes(x=fem$socM_t.avg.bright, y=fem$season.prop.ep), 
             alpha=0.5, size=3) +
  xlab("Social male throat brightness") +
  ylab("Predicted seasonal proportion of EPO") +
  ggtitle("Model predictions for female proportion EPO\nby social male throat color")

# ggsave("output-files/fem prop EPO model pred socM throat.png", h=4, w=4.5)


## social male breast-----------------------------------------------------------

# generate y-vals
pred.fem.prop.Mbreast <- custom.newdata(old.data = fem, data.type="female", scaled=T,
                                        xvar="socM_r.avg.bright_scaled", 
                                        xunscale="socM_r.avg.bright")

pred.fem.prop.Mbreast.adj1 <- calculate.glm.ci(fem.prop.adj1, pred.fem.prop.Mbreast)

ggplot(pred.fem.prop.Mbreast.adj1) +
  geom_line(aes(x=socM_r.avg.bright, y=fit), size=1.5,
            linetype="dashed") +
  geom_ribbon(aes(x=socM_r.avg.bright, ymax=upper, ymin=lower), alpha=0.3) +
  geom_point(data=fem, aes(x=fem$socM_r.avg.bright, y=fem$season.prop.ep), 
             alpha=0.5, size=3) +
  xlab("Social male breast brightness") +
  ylab("Predicted seasonal proportion of EPO") +
  ggtitle("Model predictions for female proportion EPO\nby social male throat color")

# ggsave("output-files/fem prop EPO model pred socM breast.png", h=4, w=4.5)

#-------------------------------------------------------------------------------
# Fit prop EPO for males
#-------------------------------------------------------------------------------

## Final model
male.prop.adj4 <- glm(prop.sired.epo ~ tail_scaled + socF_tail_scaled +
                        throat_avg_bright_scaled + socF_t.avg.bright_scaled +
                        breast_avg_bright_scaled + socF_r.avg.bright_scaled +
                        ci_julian_scaled, 
                      data=male, family="binomial", 
                      weights=tot.chick)

# get table of coefficients
table.male.prop.adj4 <- summary(male.prop.adj4)
table2.male.prop.adj4 <- as.data.frame(cbind(table.male.prop.adj4$coefficients, 
                              confint(male.prop.adj4)))
# back transform
table2.male.prop.adj4$beta.BT <- exp(table2.male.prop.adj4$Estimate)
table2.male.prop.adj4$lowerCI.BT <- exp(table2.male.prop.adj4$`2.5 %`)
table2.male.prop.adj4$upperCI.BT <- exp(table2.male.prop.adj4$`97.5 %`)

write.csv(table2.male.prop.adj4, 
          "output-files/male prop EPO final mod table scaled.csv")

## Mod5, add site size effect
male.prop.adj5 <- glm(prop.sired.epo ~ site.size + tail_scaled + socF_tail_scaled +
                        throat_avg_bright_scaled + socF_t.avg.bright_scaled +
                        breast_avg_bright_scaled + socF_r.avg.bright_scaled +
                        ci_julian_scaled, 
                      data=male, family="binomial", 
                      weights=tot.chick)

# get table of coefficients
table.male.prop.adj5 <- summary(male.prop.adj5)
table2.male.prop.adj5 <- as.data.frame(cbind(table.male.prop.adj5$coefficients, 
                               confint(male.prop.adj5)))
# back transform
table2.male.prop.adj5$beta.BT <- exp(table2.male.prop.adj5$Estimate)
table2.male.prop.adj5$lowerCI.BT <- exp(table2.male.prop.adj5$`2.5 %`)
table2.male.prop.adj5$upperCI.BT <- exp(table2.male.prop.adj5$`97.5 %`)

write.csv(table2.male.prop.adj5, 
          "output-files/male prop EPO final mod table site scaled.csv")

## compare models with and without site
anova(male.prop.adj4, male.prop.adj5)
# Analysis of Deviance Table
# 
# Model 1: prop.sired.epo ~ tail_scaled + socF_tail_scaled + throat_avg_bright_scaled + 
#   socF_t.avg.bright_scaled + breast_avg_bright_scaled + socF_r.avg.bright_scaled + 
#   ci_julian_scaled
# Model 2: prop.sired.epo ~ site.size + tail_scaled + socF_tail_scaled + 
#   throat_avg_bright_scaled + socF_t.avg.bright_scaled + breast_avg_bright_scaled + 
#   socF_r.avg.bright_scaled + ci_julian_scaled
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)
# 1        33     58.894                     
# 2        31     57.624  2   1.2704   0.5298


## Mod 6, full model without CI
male.prop.adj6 <- glm(prop.sired.epo ~ tail_scaled + socF_tail_scaled +
                        throat_avg_bright_scaled + socF_t.avg.bright_scaled +
                        breast_avg_bright_scaled + socF_r.avg.bright_scaled, 
                      data=male, family="binomial", 
                      weights=tot.chick)

# get table of coefficients
table.male.prop.adj6 <- summary(male.prop.adj6)
table2.male.prop.adj6 <- as.data.frame(cbind(table.male.prop.adj6$coefficients, 
                               confint(male.prop.adj6)))
# back transform
table2.male.prop.adj6$beta.BT <- exp(table2.male.prop.adj6$Estimate)
table2.male.prop.adj6$lowerCI.BT <- exp(table2.male.prop.adj6$`2.5 %`)
table2.male.prop.adj6$upperCI.BT <- exp(table2.male.prop.adj6$`97.5 %`)

write.csv(table2.male.prop.adj6, 
          "output-files/male prop EPO final table no ci scaled.csv")

###-----------------------------------------------------------------------------
# Visualize results for male propEPO
###-----------------------------------------------------------------------------

## male tail--------------------------------------------------------------------

# generate new data
pred.male.prop.tail <- custom.newdata(old.data = male, data.type="male", scaled=T,
                                      xvar="tail_scaled", 
                                      xunscale="tail")

# calculate fit and CIs
pred.male.prop.tail.adj4 <- calculate.glm.ci(male.prop.adj4, pred.male.prop.tail)

# plot
ggplot(pred.male.prop.tail.adj4) +
  geom_line(aes(x=tail, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=tail, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=tail, y=prop.sired.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Male tail length") +
  ylab("Proportion of sired offspring that are EPO")

# ggsave("output-files/male prop EPO model pred tail.png", h=4, w=4.5)


## male throat------------------------------------------------------------------

# generate new data
pred.male.prop.throat <- custom.newdata(old.data = male, data.type="male", scaled=T,
                                        xvar="throat_avg_bright_scaled", 
                                        xunscale="throat_avg_bright")

# calculate fit and CIs
pred.male.prop.throat.adj4 <- calculate.glm.ci(male.prop.adj4, 
                                               pred.male.prop.throat)

# plot
ggplot(pred.male.prop.throat.adj4) +
  geom_line(aes(x=throat_avg_bright, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=throat_avg_bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=throat_avg_bright, y=prop.sired.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Male throat average brightness") +
  ylab("Proportion of sired offspring that are EPO")

# ggsave("output-files/male prop EPO model pred throat.png", h=4, w=4.5)


## male breast -----------------------------------------------------------------

# generate y-vals
pred.male.prop.breast <- custom.newdata(old.data = male, data.type="male", scaled=T,
                                        xvar="breast_avg_bright_scaled", 
                                        xunscale="breast_avg_bright")

# calculate fit and CIs
pred.male.prop.breast.adj4 <- calculate.glm.ci(male.prop.adj4, 
                                               pred.male.prop.breast)

# plot
ggplot(pred.male.prop.breast.adj4) +
  geom_line(aes(x=breast_avg_bright, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=breast_avg_bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=breast_avg_bright, y=prop.sired.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Male breast average brightness") +
  ylab("Proportion of sired offspring that are EPO")

# ggsave("output-files/male prop EPO model pred breast.png", h=4, w=4.5)


## soc female tail--------------------------------------------------------------

# generate y-vals
pred.male.prop.Ftail <- custom.newdata(old.data = male, data.type="male", scaled=T,
                                       xvar="socF_tail_scaled", 
                                       xunscale="socF_tail")

# calculate fit and CIs
pred.male.prop.Ftail.adj4 <- calculate.glm.ci(male.prop.adj4, 
                                               pred.male.prop.Ftail)

# plot
ggplot(pred.male.prop.Ftail.adj4) +
  geom_line(aes(x=socF_tail, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=socF_tail, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=socF_tail, y=prop.sired.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social female tail length") +
  ylab("Proportion of sired offspring that are EPO")

# ggsave("output-files/male prop EPO model pred Ftail.png", h=4, w=4.5)


## soc female throat------------------------------------------------------------

# generate y-vals
pred.male.prop.Fthroat <- custom.newdata(old.data = male, data.type="male", scaled=T,
                                         xvar="socF_t.avg.bright_scaled", 
                                         xunscale="socF_t.avg.bright")

# calculate fit and CIs
pred.male.prop.Fthroat.adj4 <- calculate.glm.ci(male.prop.adj4, 
                                                pred.male.prop.Fthroat)

# plot
ggplot(pred.male.prop.Fthroat.adj4) +
  geom_line(aes(x=socF_t.avg.bright, y=fit), size=1.5, linetype="solid",
            color="darkgreen") +
  geom_ribbon(aes(x=socF_t.avg.bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=socF_t.avg.bright, y=prop.sired.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social female throat average brightness") +
  ylab("Proportion of sired offspring that are EPO")

# ggsave("output-files/male prop EPO model pred socF throat.png", h=4, w=4.5)

## soc female breast------------------------------------------------------------

# generate y-vals
pred.male.prop.Fbreast <- custom.newdata(old.data = male, data.type="male", scaled=T,
                                         xvar="socF_r.avg.bright_scaled", 
                                         xunscale="socF_r.avg.bright")
# calculate fit and CIs
pred.male.prop.Fbreast.adj4 <- calculate.glm.ci(male.prop.adj4, 
                                                pred.male.prop.Fbreast)

# plot
ggplot(pred.male.prop.Fbreast.adj4) +
  geom_line(aes(x=socF_r.avg.bright, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=socF_r.avg.bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=socF_r.avg.bright, y=prop.sired.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social female breast average brightness") +
  ylab("Proportion of sired offspring that are EPO")

# ggsave("output-files/male prop EPO model pred socF breast.png", h=4, w=4.5)

#-------------------------------------------------------------------------------
# Plot female and male results together
#-------------------------------------------------------------------------------

# set colors to create a legend
cols <- c("Females"="#CC00FF", "Males"="#00B0F0")


## Tail ------------------------------------------------------------------------

prop.focal.tail <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.prop.tail.adj1, 
            aes(x=tail, y=fit, color="Females"), 
            size=1.5, linetype="solid") +
  geom_ribbon(data=pred.fem.prop.tail.adj1, 
              aes(x=tail, ymin=lower, ymax=upper, fill="Females"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=tail, y=season.prop.ep, color="Females"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.prop.tail.adj4, 
            aes(x=tail, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.prop.tail.adj4, 
              aes(x=tail, ymax=upper, ymin=lower, fill="Males"),
              alpha=0.3) +
  geom_point(data=male, aes(x=tail, y=prop.sired.epo, color="Males"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL)+
  theme_light()

# ggsave("output-files/combined plot female and male tail on prop EPO.png", 
#        h=4, w=7, scale=0.5)



# Throat -----------------------------------------------------------------------
prop.focal.throat <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.prop.throat.adj1, 
            aes(x=throat_avg_bright, y=fit, color="Females"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.prop.throat.adj1, 
              aes(x=throat_avg_bright, ymin=lower, ymax=upper, fill="Females"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=throat_avg_bright, y=season.prop.ep, color="Females"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.prop.throat.adj4, 
            aes(x=throat_avg_bright, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.prop.throat.adj4, 
              aes(x=throat_avg_bright, ymax=upper, ymin=lower, fill="Males"),
              alpha=0.3) +
  geom_point(data=male, aes(x=throat_avg_bright, y=prop.sired.epo, color="Males"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL) +
  theme_light()

# ggsave("output-files/combined plot female and male throat on prop EPO.png", 
#        h=4, w=7, scale=0.5)



# Breast -----------------------------------------------------------------------

prop.focal.breast <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.prop.breast.adj1, 
            aes(x=breast_avg_bright, y=fit, color="Females"), 
            size=1.5, linetype="solid") +
  geom_ribbon(data=pred.fem.prop.breast.adj1, 
              aes(x=breast_avg_bright, ymin=lower, ymax=upper, fill="Females"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=breast_avg_bright, y=season.prop.ep, color="Females"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.prop.breast.adj4, 
            aes(x=breast_avg_bright, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.prop.breast.adj4, 
              aes(x=breast_avg_bright, ymax=upper, ymin=lower, fill="Males"),
              alpha=0.3) +
  geom_point(data=male, aes(x=breast_avg_bright, y=prop.sired.epo, color="Males"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL) +
 theme_light()

# ggsave("output-files/combined plot female and male breast on prop EPO.png", 
#        h=4, w=7, scale=0.5)


# Mate tail --------------------------------------------------------------------

prop.social.tail <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.prop.Mtail.adj1, 
            aes(x=socM_tail, y=fit, color="Males"), 
            size=1.5, linetype="solid") +
  geom_ribbon(data=pred.fem.prop.Mtail.adj1, 
              aes(x=socM_tail, ymin=lower, ymax=upper, fill="Males"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=socM_tail, y=season.prop.ep, color="Males"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.prop.Ftail.adj4, 
            aes(x=socF_tail, y=fit, color="Females"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.prop.Ftail.adj4, 
              aes(x=socF_tail, ymax=upper, ymin=lower, fill="Females"),
              alpha=0.3) +
  geom_point(data=male, aes(x=socF_tail, y=prop.sired.epo, color="Females"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL)+
  theme_light()

# ggsave("output-files/combined plot female and male mate tail on prop EPO.png", 
#        h=4, w=7, scale=0.5)


# Mate throat ------------------------------------------------------------------

prop.social.throat <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.prop.Mthroat.adj1, 
            aes(x=socM_t.avg.bright, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.prop.Mthroat.adj1, 
              aes(x=socM_t.avg.bright, ymin=lower, ymax=upper, fill="Males"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=socM_t.avg.bright, y=season.prop.ep, color="Males"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.prop.Fthroat.adj4, 
            aes(x=socF_t.avg.bright, y=fit, color="Females"), 
            size=1.5, linetype="solid") +
  geom_ribbon(data=pred.male.prop.Fthroat.adj4, 
              aes(x=socF_t.avg.bright, ymax=upper, ymin=lower, fill="Females"),
              alpha=0.3) +
  geom_point(data=male, aes(x=socF_t.avg.bright, y=prop.sired.epo, color="Females"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL)+
  theme_light()

# ggsave("output-files/combined plot female and male mate throat on prop EPO.png", 
#        h=4, w=7, scale=0.5)


# Mate breast ------------------------------------------------------------------

prop.social.breast <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.prop.Mbreast.adj1, 
            aes(x=socM_r.avg.bright, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.prop.Mbreast.adj1, 
              aes(x=socM_r.avg.bright, ymin=lower, ymax=upper, fill="Males"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=socM_r.avg.bright, y=season.prop.ep, color="Males"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.prop.Fbreast.adj4, 
            aes(x=socF_r.avg.bright, y=fit, color="Females"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.prop.Fbreast.adj4, 
              aes(x=socF_r.avg.bright, ymax=upper, ymin=lower, fill="Females"),
              alpha=0.3) +
  geom_point(data=male, aes(x=socF_r.avg.bright, y=prop.sired.epo, color="Females"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL)+
  theme_light()

# ggsave("output-files/combined plot female and male mate breast on prop EPO.png", 
#        h=4, w=7, scale=0.5)


## save plots to load in a different script-------------------------------------
save(prop.focal.breast, prop.focal.tail, prop.focal.throat,
     prop.social.breast, prop.social.tail, prop.social.throat,
     file="output-files/prop EP combined plots_2025-03-13.Rdata")



