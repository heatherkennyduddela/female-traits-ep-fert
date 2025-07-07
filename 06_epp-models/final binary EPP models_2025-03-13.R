
################################################################################
# Script for fitting final models for binary EPP (females and males)
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
# Fit female binomial model 
#-------------------------------------------------------------------------------

# model with raw variable units
fem.raw.mod1 <- glm(ep.yes ~ tail + socM_tail + 
                      throat_avg_bright + socM_t.avg.bright +
                      breast_avg_bright + socM_r.avg.bright +
                      ci_1_julian, data=fem, 
                    family="binomial", control=glm.control(maxit = 50))
confint(fem.raw.mod1)
exp(confint(fem.raw.mod1))*sd(fem$tail)

# final model
fem.adj.mod1 <- glm(ep.yes ~ tail_scaled + socM_tail_scaled + 
                      throat_avg_bright_scaled + socM_t.avg.bright_scaled +
                      breast_avg_bright_scaled + socM_r.avg.bright_scaled +
                      ci_1_julian_scaled, data=fem, 
                family="binomial", control=glm.control(maxit = 50))

confint(fem.adj.mod1)

# get table of coefficients
table.fem.adj.mod1 <- summary(fem.adj.mod1)
table2.fem.adj.mod1 <- as.data.frame(cbind(table.fem.adj.mod1$coefficients, 
                             confint(fem.adj.mod1)))
# Warning messages:
# 1: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# 2: glm.fit: fitted probabilities numerically 0 or 1 occurred 

# back transform
table2.fem.adj.mod1$beta.BT <- exp(table2.fem.adj.mod1$Estimate)
table2.fem.adj.mod1$lowerCI.BT <- exp(table2.fem.adj.mod1$`2.5 %`)
table2.fem.adj.mod1$upperCI.BT <- exp(table2.fem.adj.mod1$`97.5 %`)

write.csv(table2.fem.adj.mod1, "output-files/fem binary EPP final mod table scaled.csv")


# with fixed effect for site size
fem.adj.mod2 <- glm(ep.yes ~ site.size + tail_scaled + socM_tail_scaled + 
                      throat_avg_bright_scaled + socM_t.avg.bright_scaled +
                      breast_avg_bright_scaled + socM_r.avg.bright_scaled +
                      ci_1_julian_scaled, data=fem, 
                    family="binomial",
                    control=glm.control(maxit=50))

# get table of coefficients
table.fem.adj.mod2 <- summary(fem.adj.mod2)
table2.fem.adj.mod2 <- as.data.frame(cbind(table.fem.adj.mod2$coefficients, 
                             confint(fem.adj.mod2)))
# There were 50 or more warnings (use warnings() to see the first 50)
# Warning messages:
#   1: glm.fit: fitted probabilities numerically 0 or 1 occurred...

# back transform
table2.fem.adj.mod2$beta.BT <- exp(table2.fem.adj.mod2$Estimate)
table2.fem.adj.mod2$lowerCI.BT <- exp(table2.fem.adj.mod2$`2.5 %`)
table2.fem.adj.mod2$upperCI.BT <- exp(table2.fem.adj.mod2$`97.5 %`)

write.csv(table2.fem.adj.mod2, 
          "output-files/fem binary EPP final mod table site scaled.csv")

# compare models with and without site
# p-val for analysis of deviance is not significant (0.06371),
# support for dropping site
anova(fem.adj.mod1, fem.adj.mod2)

fit.both <- cbind(fem.adj.mod1$fitted.values, fem.adj.mod2$fitted.values)

## Model without CI
fem.adj.mod3 <- glm(ep.yes ~ tail_scaled + socM_tail_scaled + 
                      throat_avg_bright_scaled + socM_t.avg.bright_scaled +
                      breast_avg_bright_scaled + socM_r.avg.bright_scaled, 
                    data=fem, 
                    family="binomial", control=glm.control(maxit = 50))

# get table of coefficients
table.fem.adj.mod3 <- summary(fem.adj.mod3)
table2.fem.adj.mod3 <- as.data.frame(cbind(table.fem.adj.mod3$coefficients, 
                             confint(fem.adj.mod3)))
# back transform
table2.fem.adj.mod3$beta.BT <- exp(table2.fem.adj.mod3$Estimate)
table2.fem.adj.mod3$lowerCI.BT <- exp(table2.fem.adj.mod3$`2.5 %`)
table2.fem.adj.mod3$upperCI.BT <- exp(table2.fem.adj.mod3$`97.5 %`)

write.csv(table2.fem.adj.mod3, "output-files/fem binary EPP final table no ci scaled.csv")

###-----------------------------------------------------------------------------
#  visualize results for female Binary EP status
#-------------------------------------------------------------------------------

# NOTE: run scripts for custom functions: 
# "function for custom newdata.R"
# "function to calculate glm CIs.R"

## female tail length-----------------------------------------------------------

# generate new data
pred.fem.tail <- custom.newdata(old.data = fem, data.type="female", scaled=T,
                                xvar="tail_scaled", xunscale="tail")
# calculate fit
pred.fem.tail.adj1 <- calculate.glm.ci(fem.adj.mod1, pred.fem.tail)

ggplot(pred.fem.tail.adj1) +
  geom_line(aes(x=tail, y=fit), 
            size=1.5, color="darkgreen") +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=tail), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=tail, y=ep.yes), alpha=0.5, 
             size=3, color="darkgreen", 
             position=position_jitter(h=0.05, w=0)) +
  xlab("Female tail streamer length") +
  ylab("Probability of having EPP")

# ggsave("output-files/fem binary epp model pred fem tail.png", h=4, w=4.5)


## male tail length ------------------------------------------------------------

# new data for predictions
pred.fem.Mtail <- custom.newdata(old.data = fem, data.type="female", scaled=T,
                                 xvar="socM_tail_scaled", 
                                 xunscale="socM_tail")

# predictions and CIs from function
pred.fem.Mtail.adj1 <- calculate.glm.ci(fem.adj.mod1, pred.fem.Mtail)

# plot
ggplot(pred.fem.Mtail.adj1) +
  geom_line(aes(x=socM_tail, y=fit), size=1.5,
            color="darkgreen") +
  geom_ribbon(aes(x=socM_tail, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=fem$socM_tail, y=fem$ep.yes), alpha=0.5, size=3,
             color="darkgreen", position=position_jitter(h=0.05, w=0)) +
  xlab("Social male tail streamer length") +
  ylab("Probability of having EPP")

# ggsave("output-files/fem binary epp model pred socM tail.png", h=4, w=4.5)


## female breast color ---------------------------------------------------------

# generate y-values
pred.fem.breast <- custom.newdata(old.data = fem, data.type="female", scaled=T,
                                  xvar="breast_avg_bright_scaled", 
                                  xunscale="breast_avg_bright")

# calculate fit and CIs using function
pred.fem.breast.adj1 <- calculate.glm.ci(fem.adj.mod1, pred.fem.breast)

ggplot(pred.fem.breast.adj1) +
  geom_line(aes(x=breast_avg_bright, y=fit), size=1.5,
            color="darkgreen") +
  geom_ribbon(aes(x=breast_avg_bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=fem$breast_avg_bright, y=fem$ep.yes), 
             alpha=0.5, size=3, color="darkgreen", 
             position=position_jitter(h=0.05, w=0)) +
  xlab("Female breast brightness") +
  ylab("Probability of having EPP") 

# ggsave("output-files/fem binary epp model pred fem breast.png", h=4, w=4.5)


## Female throat color ---------------------------------------------------------

# generate new data
pred.fem.throat <- custom.newdata(old.data = fem, data.type="female", scaled=T,
                                  xvar="throat_avg_bright_scaled", 
                                  xunscale="throat_avg_bright")

# calculate fit and CIs from function
pred.fem.throat.adj1 <- calculate.glm.ci(fem.adj.mod1, pred.fem.throat)

ggplot(pred.fem.throat.adj1) +
  geom_line(aes(x=throat_avg_bright, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=throat_avg_bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=throat_avg_bright, y=ep.yes), 
             alpha=0.5, size=3, color="darkgreen", 
             position=position_jitter(h=0.05, w=0)) +
  xlab("Female throat brightness") +
  ylab("Probability of having EPP")

# ggsave("output-files/fem binary epp model pred fem throat.png", h=4, w=4.5)


# social male throat -----------------------------------------------------------

# new data for predictions
pred.fem.Mthroat <- custom.newdata(old.data = fem, data.type="female", scaled=T,
                                   xvar="socM_t.avg.bright_scaled", 
                                   xunscale="socM_t.avg.bright")

# predictions and CIs from function
pred.fem.Mthroat.adj1 <- calculate.glm.ci(fem.adj.mod1, pred.fem.Mthroat)

# plot
ggplot(pred.fem.Mthroat.adj1) +
  geom_line(aes(x=socM_t.avg.bright, y=fit), size=1.5,
            color="darkgreen", linetype="dashed") +
  geom_ribbon(aes(x=socM_t.avg.bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=fem$socM_t.avg.bright, y=fem$ep.yes), 
             alpha=0.5, size=3,
             color="darkgreen", position=position_jitter(h=0.05, w=0)) +
  xlab("Social male throat average brightness") +
  ylab("Female probability of having EPP")

# ggsave("output-files/fem binary epp model pred socM throat.png", h=4, w=4.5)


# Social male breast------------------------------------------------------------

# new data for predictions
pred.fem.Mbreast <- custom.newdata(old.data = fem, data.type="female", scaled=T,
                                   xvar="socM_r.avg.bright_scaled", 
                                   xunscale="socM_r.avg.bright")

# predictions and CIs from function
pred.fem.Mbreast.adj1 <- calculate.glm.ci(fem.adj.mod1, pred.fem.Mbreast)

# plot
ggplot(pred.fem.Mbreast.adj1) +
  geom_line(aes(x=socM_r.avg.bright, y=fit), size=1.5,
            color="darkgreen", linetype="dashed") +
  geom_ribbon(aes(x=socM_r.avg.bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=fem$socM_r.avg.bright, y=fem$ep.yes), 
             alpha=0.5, size=3,
             color="darkgreen", position=position_jitter(h=0.05, w=0)) +
  xlab("Social male breast average brightness") +
  ylab("Female probability of having EPP")

# ggsave("output-files/fem binary epp model pred socM breast.png", h=4, w=4.5)


#-------------------------------------------------------------------------------
# Fit male binomial EPP models 
#-------------------------------------------------------------------------------

# final model
male.adj.mod1 <- glm(ep.yes ~ tail_scaled + socF_tail_scaled +
                       throat_avg_bright_scaled + socF_t.avg.bright_scaled + 
                       breast_avg_bright_scaled + socF_r.avg.bright_scaled +
                       ci_julian_scaled, 
                     data=male, family="binomial",
                     control=glm.control(maxit=50))

# get table of coefficients
table.male.adj.mod1 <- summary(male.adj.mod1)
table2.male.adj.mod1 <- as.data.frame(cbind(table.male.adj.mod1$coefficients, 
                             confint(male.adj.mod1)))
# back transform
table2.male.adj.mod1$beta.BT <- exp(table2.male.adj.mod1$Estimate)
table2.male.adj.mod1$lowerCI.BT <- exp(table2.male.adj.mod1$`2.5 %`)
table2.male.adj.mod1$upperCI.BT <- exp(table2.male.adj.mod1$`97.5 %`)

write.csv(table2.male.adj.mod1, 
          "output-files/male binary EPP final mod table scaled.csv")


# with fixed effect for site size category
male.adj.mod2 <- glm(ep.yes ~ site.size + tail_scaled + socF_tail_scaled +
                       throat_avg_bright_scaled + socF_t.avg.bright_scaled + 
                       breast_avg_bright_scaled + socF_r.avg.bright_scaled +
                       ci_julian_scaled, 
                     data=male, family="binomial",
                     control=glm.control(maxit=50))

# get table of coefficients
table.male.adj.mod2 <- summary(male.adj.mod2)
table2.male.adj.mod2 <- as.data.frame(cbind(table.male.adj.mod2$coefficients, 
                              confint(male.adj.mod2)))
# back transform
table2.male.adj.mod2$beta.BT <- exp(table2.male.adj.mod2$Estimate)
table2.male.adj.mod2$lowerCI.BT <- exp(table2.male.adj.mod2$`2.5 %`)
table2.male.adj.mod2$upperCI.BT <- exp(table2.male.adj.mod2$`97.5 %`)

write.csv(table2.male.adj.mod2, 
          "output-files/male binary EPP final mod table site scaled.csv")

# compare models with and without site
# support for dropping site
anova(male.adj.mod1, male.adj.mod2)
# Analysis of Deviance Table
# 
# Model 1: ep.yes ~ tail_scaled + socF_tail_scaled + throat_avg_bright_scaled + 
#   socF_t.avg.bright_scaled + breast_avg_bright_scaled + socF_r.avg.bright_scaled + 
#   ci_julian_scaled
# Model 2: ep.yes ~ site.size + tail_scaled + socF_tail_scaled + throat_avg_bright_scaled + 
#   socF_t.avg.bright_scaled + breast_avg_bright_scaled + socF_r.avg.bright_scaled + 
#   ci_julian_scaled
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)
# 1        36     42.622                     
# 2        34     42.285  2  0.33695    0.845


# Model without CI
male.adj.mod3 <- glm(ep.yes ~ tail_scaled + socF_tail_scaled +
                       throat_avg_bright_scaled + socF_t.avg.bright_scaled + 
                       breast_avg_bright_scaled + socF_r.avg.bright_scaled, 
                     data=male, family="binomial",
                     control=glm.control(maxit=50))

# get table of coefficients
table.male.adj.mod3 <- summary(male.adj.mod3)
table2.male.adj.mod3 <- as.data.frame(cbind(table.male.adj.mod3$coefficients, 
                              confint(male.adj.mod3)))
# back transform
table2.male.adj.mod3$beta.BT <- exp(table2.male.adj.mod3$Estimate)
table2.male.adj.mod3$lowerCI.BT <- exp(table2.male.adj.mod3$`2.5 %`)
table2.male.adj.mod3$upperCI.BT <- exp(table2.male.adj.mod3$`97.5 %`)

write.csv(table2.male.adj.mod3, 
          "output-files/male binary EPP final table no ci scaled.csv")

###-----------------------------------------------------------------------------
# visualize results for male binary EP status
###-----------------------------------------------------------------------------

## female throat----------------------------------------------------------------

# generate new data
pred.male.Fthroat <- custom.newdata(old.data = male, data.type="male", scaled=T,
                                    xvar="socF_t.avg.bright_scaled", 
                                    xunscale="socF_t.avg.bright")

# calculate fit and CIs
pred.male.Fthroat.adj1 <- calculate.glm.ci(male.adj.mod1, pred.male.Fthroat)

# plot
ggplot(pred.male.Fthroat.adj1) +
  geom_line(aes(x=socF_t.avg.bright, y=fit), size=1.5, linetype="solid",
            color="darkgreen") +
  geom_ribbon(aes(x=socF_t.avg.bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) + 
  geom_point(data=male, aes(x=male$socF_t.avg.bright, y=male$ep.yes), 
             alpha=0.5, size=3, color="darkgreen", 
             position=position_jitter(h=0.05, w=0)) +
  xlab("Social female throat brightness") +
  ylab("Probability of having EPP") 

# ggsave("output-files/male binary epp model pred socF throat scaled.png", h=4, w=4.5)


## Male tail--------------------------------------------------------------------

# generate newdata
pred.male.tail <-  custom.newdata(old.data = male, data.type="male", scaled=T,
                                  xvar="tail_scaled", 
                                  xunscale="tail")

# calculate fit and CIs
pred.male.tail.adj1 <- calculate.glm.ci(male.adj.mod1, pred.male.tail)

# plot
ggplot(pred.male.tail.adj1) +
  geom_line(aes(x=tail, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=tail, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=male$tail, y=male$ep.yes), alpha=0.5, size=3,
             color="darkgreen", position=position_jitter(h=0.05, w=0)) +
  xlab("Male tail streamer length") +
  ylab("Probability of having EPP") 

# ggsave("output-files/male binary epp model pred tail.png", h=4, w=4.5)


## Male throat------------------------------------------------------------------

# generate new data
pred.male.throat <- custom.newdata(old.data = male, data.type="male", scaled=T,
                                   xvar="throat_avg_bright_scaled", 
                                   xunscale="throat_avg_bright")

# calculate fit and CIs
pred.male.throat.adj1 <- calculate.glm.ci(male.adj.mod1, pred.male.throat)

# plot
ggplot(pred.male.throat.adj1) +
  geom_line(aes(x=throat_avg_bright, y=fit), 
            size=1.5, linetype="dashed", color="darkgreen") +
  geom_ribbon(aes(x=throat_avg_bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=male$throat_avg_bright, y=male$ep.yes), 
             alpha=0.5, size=3, color="darkgreen", 
             position=position_jitter(h=0.05, w=0)) +
  xlab("Male throat brightness") +
  ylab("Probability of having EPP") 

# ggsave("output-files/male binary epp model pred throat.png", h=4, w=4.5)


## Male breast------------------------------------------------------------------

# generate new data
pred.male.breast <- custom.newdata(old.data = male, data.type="male", scaled=T,
                                   xvar="breast_avg_bright_scaled", 
                                   xunscale="breast_avg_bright")

# calculate fit and CIs
pred.male.breast.adj1 <- calculate.glm.ci(male.adj.mod1, pred.male.breast)

# plot
ggplot(pred.male.breast.adj1) +
  geom_line(aes(x=breast_avg_bright, y=fit), 
            size=1.5, linetype="dashed", color="darkgreen") +
  geom_ribbon(aes(x=breast_avg_bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=male$breast_avg_bright, y=male$ep.yes), 
             alpha=0.5, size=3, color="darkgreen",
             position=position_jitter(h=0.05, w=0)) +
  xlab("Male breast brightness") +
  ylab("Probability of having EPP")

# ggsave("output-files/male binary epp model pred breast.png", h=4, w=4.5)

# social female tail -----------------------------------------------------------

# generate new data
pred.male.Ftail <- custom.newdata(old.data = male, data.type="male", scaled=T,
                                  xvar="socF_tail_scaled", 
                                  xunscale="socF_tail")

# calculate fit and CIs
pred.male.Ftail.adj1 <- calculate.glm.ci(male.adj.mod1, pred.male.Ftail)

# plot
ggplot(pred.male.Ftail.adj1) +
  geom_line(aes(x=socF_tail, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=socF_tail, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) + 
  geom_point(data=male, aes(x=male$socF_tail, y=male$ep.yes), 
             alpha=0.5, size=3, color="darkgreen", 
             position=position_jitter(h=0.05, w=0)) +
  xlab("Social female tail streamer length") +
  ylab("Male probability of having EPP") 

# ggsave("output-files/male binary epp model pred socF tail.png", h=4, w=4.5)

# social female breast color ---------------------------------------------------

# generate new data
pred.male.Fbreast <- custom.newdata(old.data = male, data.type="male", scaled=T,
                                    xvar="socF_r.avg.bright_scaled", 
                                    xunscale="socF_r.avg.bright")

# calculate fit and CIs
pred.male.Fbreast.adj1 <- calculate.glm.ci(male.adj.mod1, pred.male.Fbreast)

# plot
ggplot(pred.male.Fbreast.adj1) +
  geom_line(aes(x=socF_r.avg.bright, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=socF_r.avg.bright, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) + 
  geom_point(data=male, aes(x=male$socF_r.avg.bright, y=male$ep.yes), 
             alpha=0.5, size=3, color="darkgreen", 
             position=position_jitter(h=0.05, w=0)) +
  xlab("Social female breast average brightness") +
  ylab("Male probability of having EPP") 

# ggsave("output-files/male binary epp model pred socF breast.png", h=4, w=4.5)


#-------------------------------------------------------------------------------
# Plot female and male results together
#-------------------------------------------------------------------------------

# set colors to create a legend
cols <- c("Females"="#F2AA84", "Males"="#0033CC")


## Tail ------------------------------------------------------------------------

binary.focal.tail <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.tail.adj1, 
            aes(x=tail, y=fit, color="Females"), 
            size=1.5, linetype="solid") +
  geom_ribbon(data=pred.fem.tail.adj1, 
               aes(x=tail, ymin=lower, ymax=upper, fill="Females"),
               alpha=0.3) +
  geom_point(data=fem, aes(x=tail, y=ep.yes, color="Females"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.tail.adj1, 
            aes(x=tail, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.tail.adj1, 
              aes(x=tail, ymax=upper, ymin=lower, fill="Males"),
              alpha=0.3) +
  geom_point(data=male, aes(x=tail, y=ep.yes, color="Males"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL) +
  theme_light()

# ggsave("output-files/combined plot female and male tail on binary EPP.png", 
#        h=4, w=7, scale=0.5)


# Throat -----------------------------------------------------------------------

binary.focal.throat <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.throat.adj1, 
            aes(x=throat_avg_bright, y=fit, color="Females"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.throat.adj1, 
              aes(x=throat_avg_bright, ymin=lower, ymax=upper, fill="Females"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=throat_avg_bright, y=ep.yes, color="Females"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.throat.adj1, 
            aes(x=throat_avg_bright, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.throat.adj1, 
              aes(x=throat_avg_bright, ymax=upper, ymin=lower, fill="Males"),
              alpha=0.3) +
  geom_point(data=male, aes(x=throat_avg_bright, y=ep.yes, color="Males"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL) +
  theme_light()


# ggsave("output-files/combined plot female and male throat on binary EPP.png", 
#        h=4, w=7, scale=0.5)

# Breast -----------------------------------------------------------------------

binary.focal.breast <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.breast.adj1, 
            aes(x=breast_avg_bright, y=fit, color="Females"), 
            size=1.5, linetype="solid") +
  geom_ribbon(data=pred.fem.breast.adj1, 
              aes(x=breast_avg_bright, ymin=lower, ymax=upper, fill="Females"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=breast_avg_bright, y=ep.yes, color="Females"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.breast.adj1, 
            aes(x=breast_avg_bright, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.breast.adj1, 
              aes(x=breast_avg_bright, ymax=upper, ymin=lower, fill="Males"),
              alpha=0.3) +
  geom_point(data=male, aes(x=breast_avg_bright, y=ep.yes, color="Males"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL) +
  theme_light()


# ggsave("output-files/combined plot female and male breast on binary EPP.png", 
#        h=4, w=7, scale=0.5)


## social mate tail-------------------------------------------------------------

binary.social.tail <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.Mtail.adj1, 
            aes(x=socM_tail, y=fit, color="Males"), 
            size=1.5, linetype="solid") +
  geom_ribbon(data=pred.fem.Mtail.adj1, 
              aes(x=socM_tail, ymin=lower, ymax=upper, fill="Males"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=socM_tail, y=ep.yes, color="Males"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.Ftail.adj1, 
            aes(x=socF_tail, y=fit, color="Females"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.Ftail.adj1, 
              aes(x=socF_tail, ymax=upper, ymin=lower, fill="Females"),
              alpha=0.3) +
  geom_point(data=male, aes(x=socF_tail, y=ep.yes, color="Females"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL)+
  theme_light()

# ggsave("output-files/combined plot female and male mate tail on binary EPP.png", 
#        h=4, w=7, scale=0.5)


## social mate throat-----------------------------------------------------------

binary.social.throat <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.Mthroat.adj1, 
            aes(x=socM_t.avg.bright, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.Mthroat.adj1, 
              aes(x=socM_t.avg.bright, ymin=lower, ymax=upper, fill="Males"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=socM_t.avg.bright, y=ep.yes, color="Males"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.Fthroat.adj1, 
            aes(x=socF_t.avg.bright, y=fit, color="Females"), 
            size=1.5, linetype="solid") +
  geom_ribbon(data=pred.male.Fthroat.adj1, 
              aes(x=socF_t.avg.bright, ymax=upper, ymin=lower, fill="Females"),
              alpha=0.3) +
  geom_point(data=male, aes(x=socF_t.avg.bright, y=ep.yes, color="Females"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL)+
  theme_light()

# ggsave("output-files/combined plot female and male mate throat on binary EPP.png", 
#        h=4, w=7, scale=0.5)


## social mate breast-----------------------------------------------------------

binary.social.breast <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.Mbreast.adj1, 
            aes(x=socM_r.avg.bright, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.Mbreast.adj1, 
              aes(x=socM_r.avg.bright, ymin=lower, ymax=upper, fill="Males"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=socM_r.avg.bright, y=ep.yes, color="Males"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.Fbreast.adj1, 
            aes(x=socF_r.avg.bright, y=fit, color="Females"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.Fbreast.adj1, 
              aes(x=socF_r.avg.bright, ymax=upper, ymin=lower, fill="Females"),
              alpha=0.3) +
  geom_point(data=male, aes(x=socF_r.avg.bright, y=ep.yes, color="Females"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL)+
  theme_light()

# ggsave("output-files/combined plot female and male mate breast on binary EPP.png", 
#        h=4, w=7, scale=0.5)


## save plots to load in a different script-------------------------------------
save(binary.focal.breast, binary.focal.tail, binary.focal.throat, 
     binary.social.breast, binary.social.tail, binary.social.throat,
     file="output-files/binary EP combined plots_2025-03-13_bluetan.Rdata")




