
################################################################################
# Script for fitting models of number of mates for females and males
# Heather Kenny-Duddela
# Feb 6, 2025
################################################################################

# load data
fem.season <- read.csv("input-files/female_season_table_2025-03-13.csv")

male.season <- read.csv("input-files/male_season_table_2025-03-13.csv")

# libraries
library(tidyverse) #ggplot2 and dplyr
library(car) # for Anova and vif (variance inflation factor)
library(broom)# for calculating confidence intervals around coefficient estimates
library(AICcmodavg) # for calculating AIC
library(MASS) # for the glm.nb function

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
# Fit linear models for females
#-------------------------------------------------------------------------------

fem.numMates.adj1 <- lm(num.tot.mates ~ tail_scaled + socM_tail_scaled + 
                            throat_avg_bright_scaled + socM_t.avg.bright_scaled +
                            breast_avg_bright_scaled + socM_r.avg.bright_scaled +
                            ci_1_julian_scaled, data=fem)

# with fixed effect of site size category
fem.numMates.adj2 <- lm(num.tot.mates ~ site.size + tail_scaled + socM_tail_scaled + 
                            throat_avg_bright_scaled + socM_t.avg.bright_scaled +
                            breast_avg_bright_scaled + socM_r.avg.bright_scaled +
                            ci_1_julian_scaled, data=fem)

# get table of coefficients
table.fem.numMates.adj1 <- summary(fem.numMates.adj1)
table2.fem.numMates.adj1 <- as.data.frame(cbind(table.fem.numMates.adj1$coefficients, 
                                confint(fem.numMates.adj1)))

write.csv(table2.fem.numMates.adj1, "output-files/fem numMates final mod table scale.csv")

# get table of coefficients
table.fem.numMates.adj2 <- summary(fem.numMates.adj2)
table2.fem.numMates.adj2 <- as.data.frame(cbind(table.fem.numMates.adj2$coefficients, 
                                confint(fem.numMates.adj2)))
# back transform
table2.fem.numMates.adj2$beta.BT <- exp(table2.fem.numMates.adj2$Estimate)
table2.fem.numMates.adj2$lowerCI.BT <- exp(table2.fem.numMates.adj2$`2.5 %`)
table2.fem.numMates.adj2$upperCI.BT <- exp(table2.fem.numMates.adj2$`97.5 %`)

write.csv(table2.fem.numMates.adj2, "output-files/fem numMates final mod table site scale.csv")

## compare models with and without site
# support for dropping site
anova(fem.numMates.adj1, fem.numMates.adj2)
# Analysis of Variance Table
# 
# Model 1: num.tot.mates ~ tail + socM_tail + throat_avg_bright + socM_t.avg.bright + 
#   breast_avg_bright + socM_r.avg.bright + ci_1_julian
# Model 2: num.tot.mates ~ site.size + tail + socM_tail + throat_avg_bright + 
#   socM_t.avg.bright + breast_avg_bright + socM_r.avg.bright + 
#   ci_1_julian
# Res.Df    RSS Df Sum of Sq      F Pr(>F)
# 1     35 25.095                           
# 2     33 24.898  2   0.19636 0.1301 0.8784


## Model without CI
fem.numMates.adj3 <- lm(num.tot.mates ~ tail_scaled + socM_tail_scaled + 
                            throat_avg_bright_scaled + socM_t.avg.bright_scaled +
                            breast_avg_bright_scaled + socM_r.avg.bright_scaled,
                          data=fem)
# get table of coefficients
table.fem.numMates.adj3 <- summary(fem.numMates.adj3)
table2.fem.numMates.adj3 <- as.data.frame(cbind(table.fem.numMates.adj3$coefficients, 
                                confint(fem.numMates.adj3)))
# back transform
table2.fem.numMates.adj3$beta.BT <- exp(table2.fem.numMates.adj3$Estimate)
table2.fem.numMates.adj3$lowerCI.BT <- exp(table2.fem.numMates.adj3$`2.5 %`)
table2.fem.numMates.adj3$upperCI.BT <- exp(table2.fem.numMates.adj3$`97.5 %`)

write.csv(table2.fem.numMates.adj3, "output-files/fem numMates final table no ci scaled.csv")

#-------------------------------------------------------------------------------
# Visualize results for female model adj1
#-------------------------------------------------------------------------------

# NOTE: run scripts for custom functions: 
# "function for custom newdata.R"
# "function to calculate glm CIs.R"

## female tail length-----------------------------------------------------------

# generate new data
pred.fem.numMates.tail <- custom.newdata(old.data = fem, data.type="female", scaled=T,
                                         xvar="tail_scaled", xunscale="tail")
# calculate fit
results.tail.adj1 <- predict.lm(fem.numMates.adj1, newdata = pred.fem.numMates.tail, 
                      interval="confidence", level=0.95)
# combine tables
pred.fem.numMates.tail.adj1 <- cbind(pred.fem.numMates.tail, results.tail.adj1)

# plot
ggplot(pred.fem.numMates.tail.adj1) +
  geom_line(aes(x=tail, y=fit), linewidth=1.5, linetype="solid", color="darkgreen") +
  geom_ribbon(aes(x=tail, ymax=upr, ymin=lwr), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=tail, y=num.tot.mates), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Female tail length") +
  ylab("Number of genetic mates")


## female throat color----------------------------------------------------------

# vector for range of throat color
range(fem$throat, na.rm=T) # 9,25
seq.fem.throat <- seq(9, 25, 0.01)

# generate y-values
pred.fem.numMates.throat <- custom.newdata(old.data = fem, data.type="female", 
                                           scaled=T,
                                           xvar="throat_avg_bright_scaled", 
                                           xunscale="throat_avg_bright")

results.throat.adj1 <- predict.lm(fem.numMates.adj1, newdata = pred.fem.numMates.throat, 
                                interval="confidence", level=0.95)

pred.fem.numMates.throat.adj1 <- cbind(pred.fem.numMates.throat, results.throat.adj1)

ggplot(pred.fem.numMates.throat.adj1) +
  geom_line(aes(x=throat_avg_bright, y=fit), linewidth=1.5, linetype="dashed", color="darkgreen") +
  geom_ribbon(aes(x=throat_avg_bright, ymax=upr, ymin=lwr), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=throat_avg_bright, y=num.tot.mates), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Female throat brightness") +
  ylab("Number of genetic mates")


## female breast color----------------------------------------------------------

# generate new data
pred.fem.numMates.breast <- custom.newdata(old.data = fem, data.type="female", 
                                           scaled=T,
                                           xvar="breast_avg_bright_scaled", 
                                           xunscale="breast_avg_bright")

results.breast.adj1 <- predict.lm(fem.numMates.adj1, newdata = pred.fem.numMates.breast, 
                                  interval="confidence", level=0.95)

pred.fem.numMates.breast.adj1 <- cbind(pred.fem.numMates.breast, results.breast.adj1)

ggplot(pred.fem.numMates.breast.adj1) +
  geom_line(aes(x=breast_avg_bright, y=fit), linewidth=1.5, linetype="dashed", color="darkgreen") +
  geom_ribbon(aes(x=breast_avg_bright, ymax=upr, ymin=lwr), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=breast_avg_bright, y=num.tot.mates), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Female breast brightness") +
  ylab("Number of genetic mates")


## social male tail length------------------------------------------------------

# generate new data
pred.fem.numMates.Mtail <- custom.newdata(old.data = fem, data.type="female", 
                                          scaled=T,
                                          xvar="socM_tail_scaled", 
                                          xunscale="socM_tail")

results.Mtail.adj1 <- predict.lm(fem.numMates.adj1, newdata = pred.fem.numMates.Mtail, 
                                  interval="confidence", level=0.95)

pred.fem.numMates.Mtail.adj1 <- cbind(pred.fem.numMates.Mtail, results.Mtail.adj1)

ggplot(pred.fem.numMates.Mtail.adj1) +
  geom_line(aes(x=socM_tail, y=fit), linewidth=1.5, linetype="dashed", color="darkgreen") +
  geom_ribbon(aes(x=socM_tail, ymax=upr, ymin=lwr), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=socM_tail, y=num.tot.mates), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social male tail length") +
  ylab("Number of genetic mates")


## social male throat-----------------------------------------------------------

# generate y-values
pred.fem.numMates.Mthroat <- custom.newdata(old.data = fem, data.type="female", 
                                            scaled=T,
                                            xvar="socM_t.avg.bright_scaled", 
                                            xunscale="socM_t.avg.bright")

results.Mthroat.adj1 <- predict.lm(fem.numMates.adj1, newdata = pred.fem.numMates.Mthroat, 
                                 interval="confidence", level=0.95)

pred.fem.numMates.Mthroat.adj1 <- cbind(pred.fem.numMates.Mthroat, results.Mthroat.adj1)

ggplot(pred.fem.numMates.Mthroat.adj1) +
  geom_line(aes(x=socM_t.avg.bright, y=fit), linewidth=1.5, linetype="dashed", color="darkgreen") +
  geom_ribbon(aes(x=socM_t.avg.bright, ymax=upr, ymin=lwr), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=socM_t.avg.bright, y=num.tot.mates), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social male throat brightness") +
  ylab("Number of genetic mates")


## social male breast color-----------------------------------------------------

# generate y-values
pred.fem.numMates.Mbreast <- custom.newdata(old.data = fem, data.type="female", 
                                             scaled=T,
                                             xvar="socM_r.avg.bright_scaled", 
                                             xunscale="socM_r.avg.bright")

results.Mbreast.adj1 <- predict.lm(fem.numMates.adj1, newdata = pred.fem.numMates.Mbreast, 
                                   interval="confidence", level=0.95)

pred.fem.numMates.Mbreast.adj1 <- cbind(pred.fem.numMates.Mbreast, results.Mbreast.adj1)

ggplot(pred.fem.numMates.Mbreast.adj1) +
  geom_line(aes(x=socM_r.avg.bright, y=fit), linewidth=1.5, linetype="dashed", color="darkgreen") +
  geom_ribbon(aes(x=socM_r.avg.bright, ymax=upr, ymin=lwr), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=socM_r.avg.bright, y=num.tot.mates), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social male breast brightness") +
  ylab("Number of genetic mates")


################################################################################

#-------------------------------------------------------------------------------
# Fit Poisson model for males
#-------------------------------------------------------------------------------

male.numMates.adj1 <- glm(num.tot.mates ~ tail_scaled + socF_tail_scaled +
                             throat_avg_bright_scaled + socF_t.avg.bright_scaled + 
                             breast_avg_bright_scaled + socF_r.avg.bright_scaled +
                             ci_julian_scaled, data=male, family="poisson")

# with site size effect
male.numMates.adj2 <- glm(num.tot.mates ~ site.size + tail_scaled + socF_tail_scaled +
                             throat_avg_bright_scaled + socF_t.avg.bright_scaled + 
                             breast_avg_bright_scaled + socF_r.avg.bright_scaled +
                             ci_julian_scaled, data=male,family="poisson")

# compare models with and without site
# support for dropping site
anova(male.numMates.adj1, male.numMates.adj2)
# Analysis of Deviance Table
# 
# Model 1: num.tot.mates ~ tail_scaled + socF_tail_scaled + throat_avg_bright_scaled + 
#   socF_t.avg.bright_scaled + breast_avg_bright_scaled + socF_r.avg.bright_scaled + 
#   ci_julian_scaled
# Model 2: num.tot.mates ~ site.size + tail_scaled + socF_tail_scaled + 
#   throat_avg_bright_scaled + socF_t.avg.bright_scaled + breast_avg_bright_scaled + 
#   socF_r.avg.bright_scaled + ci_julian_scaled
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)
# 1        36     21.803                     
# 2        34     20.471  2   1.3313   0.5139

# get table of coefficients
table.male.numMates.adj1 <- summary(male.numMates.adj1)
table2.male.numMates.adj1 <- as.data.frame(cbind(table.male.numMates.adj1$coefficients, 
                                 confint(male.numMates.adj1)))
# back transform
table2.male.numMates.adj1$beta.BT <- exp(table2.male.numMates.adj1$Estimate)
table2.male.numMates.adj1$lowerCI.BT <- exp(table2.male.numMates.adj1$`2.5 %`)
table2.male.numMates.adj1$upperCI.BT <- exp(table2.male.numMates.adj1$`97.5 %`)

write.csv(table2.male.numMates.adj1, "output-files/male numMates final mod table scaled.csv")

# get table of coefficients
table.male.numMates.adj2 <- summary(male.numMates.adj2)
table2.male.numMates.adj2 <- as.data.frame(cbind(table.male.numMates.adj2$coefficients, 
                                 confint(male.numMates.adj2)))
# back transform
table2.male.numMates.adj2$beta.BT <- exp(table2.male.numMates.adj2$Estimate)
table2.male.numMates.adj2$lowerCI.BT <- exp(table2.male.numMates.adj2$`2.5 %`)
table2.male.numMates.adj2$upperCI.BT <- exp(table2.male.numMates.adj2$`97.5 %`)

write.csv(table2.male.numMates.adj2, "output-files/male numMates final mod table site scaled.csv")


## Model without CI
male.numMates.adj3 <- glm(num.tot.mates ~ tail_scaled + socF_tail_scaled +
                             throat_avg_bright_scaled + socF_t.avg.bright_scaled + 
                             breast_avg_bright_scaled + socF_r.avg.bright_scaled,
                           data=male, family="poisson")
# get table of coefficients
table.male.numMates.adj3 <- summary(male.numMates.adj3)
table2.male.numMates.adj3 <- as.data.frame(cbind(table.male.numMates.adj3$coefficients, 
                                 confint(male.numMates.adj3)))
# back transform
table2.male.numMates.adj3$beta.BT <- exp(table2.male.numMates.adj3$Estimate)
table2.male.numMates.adj3$lowerCI.BT <- exp(table2.male.numMates.adj3$`2.5 %`)
table2.male.numMates.adj3$upperCI.BT <- exp(table2.male.numMates.adj3$`97.5 %`)

write.csv(table2.male.numMates.adj3, "output-files/male numMates final table no ci scaled.csv")

#-------------------------------------------------------------------------------
# Visualize results for male model adj1
#-------------------------------------------------------------------------------

## social female throat color---------------------------------------------------

# generate new data
pred.male.numMates.Fthroat <- custom.newdata(old.data = male, data.type="male", 
                                             scaled=T,
                                             xvar="socF_t.avg.bright_scaled", 
                                             xunscale="socF_t.avg.bright")

# calculate fit and CIs
pred.male.numMates.Fthroat.adj1 <- calculate.glm.ci(male.numMates.adj1, 
                                                  pred.male.numMates.Fthroat)

ggplot(pred.male.numMates.Fthroat.adj1) +
  geom_line(aes(x=socF_t.avg.bright, y=fit), linewidth=1.5, linetype="solid",
            color="darkgreen") +
  geom_ribbon(aes(x=socF_t.avg.bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=socF_t.avg.bright, y=num.tot.mates), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social female throat average brightness") +
  ylab("Male number of genetic dams")


## social female tail-----------------------------------------------------------

# generate new data
pred.male.numMates.Ftail <- custom.newdata(old.data = male, data.type="male", 
                                           scaled=T,
                                           xvar="socF_tail_scaled", 
                                           xunscale="socF_tail")

# calculate fit and CIs
pred.male.numMates.Ftail.adj1 <- calculate.glm.ci(male.numMates.adj1, 
                                                pred.male.numMates.Ftail)

ggplot(pred.male.numMates.Ftail.adj1) +
  geom_line(aes(x=socF_tail, y=fit), linewidth=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=socF_tail, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=socF_tail, y=num.tot.mates), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social female tail streamer length") +
  ylab("Male number of genetic dams")


## Social female breast---------------------------------------------------------

# generate new data
pred.male.numMates.Fbreast <- custom.newdata(old.data = male, data.type="male", 
                                             scaled=T,
                                             xvar="socF_r.avg.bright_scaled", 
                                             xunscale="socF_r.avg.bright")

# calculate fit and CIs
pred.male.numMates.Fbreast.adj1 <- calculate.glm.ci(male.numMates.adj1, 
                                                  pred.male.numMates.Fbreast)

ggplot(pred.male.numMates.Fbreast.adj1) +
  geom_line(aes(x=socF_r.avg.bright, y=fit), linewidth=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=socF_r.avg.bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=socF_r.avg.bright, y=num.tot.mates), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social female breast average brightness") +
  ylab("Male number of genetic dams")

## male tail--------------------------------------------------------------------

# generate new data
pred.male.numMates.tail <- custom.newdata(old.data = male, data.type="male", 
                                          scaled=T,
                                          xvar="tail_scaled", 
                                          xunscale="tail")
# calculate CI and fit
pred.male.numMates.tail.adj1 <- calculate.glm.ci(male.numMates.adj1,
                                               pred.male.numMates.tail)
ggplot(pred.male.numMates.tail.adj1) +
  geom_line(aes(x=tail, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=tail, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=tail, y=num.tot.mates), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Male tail length") +
  ylab("Number of genetic dams")


## Male throat------------------------------------------------------------------ 

# generate y-values
pred.male.numMates.throat <- custom.newdata(old.data = male, data.type="male", 
                                            scaled=T,
                                            xvar="throat_avg_bright_scaled", 
                                            xunscale="throat_avg_bright")

# calculate fit and CIs
pred.male.numMates.throat.adj1 <- calculate.glm.ci(male.numMates.adj1,
                                                 pred.male.numMates.throat)
ggplot(pred.male.numMates.throat.adj1) +
  geom_line(aes(x=throat_avg_bright, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=throat_avg_bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=throat_avg_bright, y=num.tot.mates), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Male throat average brightness") +
  ylab("Male number of genetic dams")


## Male breast------------------------------------------------------------------

# generate y-values
pred.male.numMates.breast <- custom.newdata(old.data = male, data.type="male", 
                                            scaled=T,
                                            xvar="breast_avg_bright_scaled", 
                                            xunscale="breast_avg_bright")

# calculate fit and CIs
pred.male.numMates.breast.adj1 <- calculate.glm.ci(male.numMates.adj1, 
                                                 pred.male.numMates.breast)

ggplot(pred.male.numMates.breast.adj1) +
  geom_line(aes(x=breast_avg_bright, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=breast_avg_bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=breast_avg_bright, y=num.tot.mates), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Male breast average brightness") +
  ylab("Male number of genetic dams")


#-------------------------------------------------------------------------------
# # Visualize female and male numMates together
#-------------------------------------------------------------------------------

# set colors to create a legend
cols <- c("Females"="#F2AA84", "Males"="#0033CC")


## Tail ------------------------------------------------------------------------

focal.tail <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.numMates.tail.adj1, 
            aes(x=tail, y=fit, color="Females"), 
            size=1.5, linetype="solid") +
  geom_ribbon(data=pred.fem.numMates.tail.adj1, 
              aes(x=tail, ymin=lwr, ymax=upr, fill="Females"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=tail, y=num.tot.mates, color="Females"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.numMates.tail.adj1, 
            aes(x=tail, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.numMates.tail.adj1, 
              aes(x=tail, ymax=upper, ymin=lower, fill="Males"),
              alpha=0.3) +
  geom_point(data=male, aes(x=tail, y=num.tot.mates, color="Males"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL) +
  theme_light()

# ggsave("output-files/combined plot female and male tail on numMates.png", 
#        h=4, w=7, scale=0.5)


# Throat -----------------------------------------------------------------------

focal.throat <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.numMates.throat.adj1, 
            aes(x=throat_avg_bright, y=fit, color="Females"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.numMates.throat.adj1, 
              aes(x=throat_avg_bright, ymin=lwr, ymax=upr, fill="Females"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=throat_avg_bright, y=num.tot.mates, color="Females"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.numMates.throat.adj1, 
            aes(x=throat_avg_bright, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.numMates.throat.adj1, 
              aes(x=throat_avg_bright, ymax=upper, ymin=lower, fill="Males"),
              alpha=0.3) +
  geom_point(data=male, aes(x=throat_avg_bright, y=num.tot.mates, color="Males"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL) +
  theme_light()

# ggsave("output-files/combined plot female and male throat on numMates.png", 
#        h=4, w=7, scale=0.5)

# Breast -----------------------------------------------------------------------

focal.breast <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.numMates.breast.adj1, 
            aes(x=breast_avg_bright, y=fit, color="Females"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.numMates.breast.adj1, 
              aes(x=breast_avg_bright, ymin=lwr, ymax=upr, fill="Females"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=breast_avg_bright, y=num.tot.mates, color="Females"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.numMates.breast.adj1, 
            aes(x=breast_avg_bright, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.numMates.breast.adj1, 
              aes(x=breast_avg_bright, ymax=upper, ymin=lower, fill="Males"),
              alpha=0.3) +
  geom_point(data=male, aes(x=breast_avg_bright, y=num.tot.mates, color="Males"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL) +
  theme_light()

# ggsave("output-files/combined plot female and male breast on numMates.png", 
#       h=4, w=7, scale=0.5)



## social mate tail-------------------------------------------------------------

social.tail <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.numMates.Mtail.adj1, 
            aes(x=socM_tail, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.numMates.Mtail.adj1, 
              aes(x=socM_tail, ymin=lwr, ymax=upr, fill="Males"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=socM_tail, y=num.tot.mates, color="Males"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.numMates.Ftail.adj1, 
            aes(x=socF_tail, y=fit, color="Females"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.numMates.Ftail.adj1, 
              aes(x=socF_tail, ymax=upper, ymin=lower, fill="Females"),
              alpha=0.3) +
  geom_point(data=male, aes(x=socF_tail, y=num.tot.mates, color="Females"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL)+
  theme_light()

# ggsave("output-files/combined plot female and male mate tail on numMates.png", 
#        h=4, w=7, scale=0.5)


## social mate throat-----------------------------------------------------------

social.throat <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.numMates.Mthroat.adj1, 
            aes(x=socM_t.avg.bright, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.numMates.Mthroat.adj1, 
              aes(x=socM_t.avg.bright, ymin=lwr, ymax=upr, fill="Males"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=socM_t.avg.bright, y=num.tot.mates, color="Males"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.numMates.Fthroat.adj1, 
            aes(x=socF_t.avg.bright, y=fit, color="Females"), 
            size=1.5, linetype="solid") +
  geom_ribbon(data=pred.male.numMates.Fthroat.adj1, 
              aes(x=socF_t.avg.bright, ymax=upper, ymin=lower, fill="Females"),
              alpha=0.3) +
  geom_point(data=male, aes(x=socF_t.avg.bright, y=num.tot.mates, color="Females"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL)+
  theme_light()

# ggsave("output-files/combined plot female and male mate throat on numMates.png", 
#        h=4, w=7, scale=0.5)


## social mate breast-----------------------------------------------------------

social.breast <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.numMates.Mbreast.adj1, 
            aes(x=socM_r.avg.bright, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.numMates.Mbreast.adj1, 
              aes(x=socM_r.avg.bright, ymin=lwr, ymax=upr, fill="Males"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=socM_r.avg.bright, y=num.tot.mates, color="Males"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.numMates.Fbreast.adj1, 
            aes(x=socF_r.avg.bright, y=fit, color="Females"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.numMates.Fbreast.adj1, 
              aes(x=socF_r.avg.bright, ymax=upper, ymin=lower, fill="Females"),
              alpha=0.3) +
  geom_point(data=male, aes(x=socF_r.avg.bright, y=num.tot.mates, color="Females"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL)+
  theme_light()

# ggsave("output-files/combined plot female and male mate breast on numMates.png", 
#        h=4, w=7, scale=0.5)


## save plots as Rdata to load in a different script----------------------------

save(focal.tail, focal.throat, focal.breast, social.tail, social.breast,
     social.throat, file="output-files/numMates combined plots_2025-03-13_bluetan.Rdata")


