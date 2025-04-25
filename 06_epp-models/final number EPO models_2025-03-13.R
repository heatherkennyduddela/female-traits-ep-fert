
################################################################################
# Script for fitting final models for number of EPO
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
# Fit female model number EPO
#-------------------------------------------------------------------------------

fem.numEPO.adj1 <- glm.nb(num.epo ~ tail_scaled + socM_tail_scaled + 
                            throat_avg_bright_scaled + socM_t.avg.bright_scaled +
                            breast_avg_bright_scaled + socM_r.avg.bright_scaled +
                            ci_1_julian_scaled, data=fem)

# with fixed effect of site size category
fem.numEPO.adj2 <- glm.nb(num.epo ~ site.size + tail_scaled + socM_tail_scaled + 
                            throat_avg_bright_scaled + socM_t.avg.bright_scaled +
                            breast_avg_bright_scaled + socM_r.avg.bright_scaled +
                            ci_1_julian_scaled, data=fem,
                          control=glm.control(maxit = 50))
summary(fem.numEPO.adj2)

# get table of coefficients
table.fem.numEPO.adj1 <- summary(fem.numEPO.adj1)
table2.fem.numEPO.adj1 <- as.data.frame(cbind(table.fem.numEPO.adj1$coefficients, 
                             confint(fem.numEPO.adj1)))
# back-transform
table2.fem.numEPO.adj1$beta.BT <- exp(table2.fem.numEPO.adj1$Estimate)
table2.fem.numEPO.adj1$lowerCI.BT <- exp(table2.fem.numEPO.adj1$`2.5 %`)
table2.fem.numEPO.adj1$upperCI.BT <- exp(table2.fem.numEPO.adj1$`97.5 %`)

write.csv(table2.fem.numEPO.adj1, "output-files/fem numEPO final mod table scaled.csv")

# get table of coefficients
table.fem.numEPO.adj2 <- summary(fem.numEPO.adj2)
table2.fem.numEPO.adj2 <- as.data.frame(cbind(table.fem.numEPO.adj2$coefficients, 
                                confint(fem.numEPO.adj2)))
# back-transform
table2.fem.numEPO.adj2$beta.BT <- exp(table2.fem.numEPO.adj2$Estimate)
table2.fem.numEPO.adj2$lowerCI.BT <- exp(table2.fem.numEPO.adj2$`2.5 %`)
table2.fem.numEPO.adj2$upperCI.BT <- exp(table2.fem.numEPO.adj2$`97.5 %`)

write.csv(table2.fem.numEPO.adj2, "output-files/fem numEPO final mod table site scaled.csv")

## compare models with and without site
# support for dropping site
anova(fem.numEPO.adj1, fem.numEPO.adj2)

## Model without CI
fem.numEPO.adj3 <- glm.nb(num.epo ~ tail_scaled + socM_tail_scaled + 
                            throat_avg_bright_scaled + socM_t.avg.bright_scaled +
                            breast_avg_bright_scaled + socM_r.avg.bright_scaled,
                          data=fem)
# get table of coefficients
table.fem.numEPO.adj3 <- summary(fem.numEPO.adj3)
table2.fem.numEPO.adj3 <- as.data.frame(cbind(table.fem.numEPO.adj3$coefficients, 
                                confint(fem.numEPO.adj3)))
# back-transform
table2.fem.numEPO.adj3$beta.BT <- exp(table2.fem.numEPO.adj3$Estimate)
table2.fem.numEPO.adj3$lowerCI.BT <- exp(table2.fem.numEPO.adj3$`2.5 %`)
table2.fem.numEPO.adj3$upperCI.BT <- exp(table2.fem.numEPO.adj3$`97.5 %`)

write.csv(table2.fem.numEPO.adj3, "output-files/fem numEPO final table no ci scaled.csv")

#-------------------------------------------------------------------------------
###  visualize results for adj1 model
#-------------------------------------------------------------------------------

## female tail length-----------------------------------------------------------

# vector for range of tail lengths
seq.fem.tail <- seq(range(fem$tail_scaled)[1], range(fem$tail_scaled)[2], 0.005)
fem.tail.n <- length(seq.fem.tail)

# generate new data
pred.fem.numEPO.tail <- 
  data.frame(tail_scaled = seq.fem.tail,
             site.size = rep("medium", fem.tail.n),
       socM_tail_scaled = rep(mean(fem$socM_tail_scaled, na.rm=T), fem.tail.n),
       ci_1_julian_scaled = rep(mean(fem$ci_1_julian_scaled, na.rm=T), fem.tail.n),
       throat_avg_bright_scaled = rep(mean(fem$throat_avg_bright_scaled, na.rm=T), fem.tail.n),
       socM_t.avg.bright_scaled = rep(mean(fem$socM_t.avg.bright_scaled, na.rm=T), fem.tail.n),
       breast_avg_bright_scaled = rep(mean(fem$breast_avg_bright_scaled, na.rm=T), fem.tail.n),
       socM_r.avg.bright_scaled = rep(mean(fem$socM_r.avg.bright_scaled, na.rm=T), fem.tail.n))

# calculate fit and CIs
pred.fem.numEPO.tail.adj1 <- calculate.glm.ci(fem.numEPO.adj1, 
                                              pred.fem.numEPO.tail)
# un-scale the predictor
pred.fem.numEPO.tail.adj1$tail.unscale <- (pred.fem.numEPO.tail.adj1$tail_scaled *
                                             sd(fem$tail) + mean(fem$tail))

ggplot(pred.fem.numEPO.tail.adj1) +
  geom_line(aes(x=tail.unscale, y=fit), size=1.5, linetype="solid", color="darkgreen") +
  geom_ribbon(aes(x=tail.unscale, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=tail, y=num.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Female tail length") +
  ylab("Number of EPO")

# ggsave("output-files/fem num EPO model pred fem tail.png", h=4, w=4.5)


## Female throat color----------------------------------------------------------

# vector for range of throat color
seq.fem.throat <- seq(range(fem$throat_avg_bright_scaled)[1], 
                      range(fem$throat_avg_bright_scaled)[2], 0.001)
fem.throat.n <- length(seq.fem.throat)

# generate y-values
pred.fem.numEPO.throat <-  
  data.frame(tail_scaled = rep(mean(fem$tail_scaled, narm=T), fem.throat.n),
             site.size = rep("medium", fem.throat.n),
       socM_tail_scaled = rep(mean(fem$socM_tail_scaled, na.rm=T), fem.throat.n),
       ci_1_julian_scaled = rep(mean(fem$ci_1_julian_scaled, na.rm=T), fem.throat.n),
       throat_avg_bright_scaled = seq.fem.throat,
       socM_t.avg.bright_scaled = rep(mean(fem$socM_t.avg.bright_scaled, na.rm=T), fem.throat.n),
       breast_avg_bright_scaled = rep(mean(fem$breast_avg_bright_scaled, na.rm=T), fem.throat.n),
       socM_r.avg.bright_scaled = rep(mean(fem$socM_r.avg.bright_scaled, na.rm=T), fem.throat.n))

# calculate fit and CIs
pred.fem.numEPO.throat.adj1 <- calculate.glm.ci(fem.numEPO.adj1, 
                                                pred.fem.numEPO.throat)
# un-scale the predictor
pred.fem.numEPO.throat.adj1$throat.unscale <- (pred.fem.numEPO.throat.adj1$throat_avg_bright_scaled *
                                             sd(fem$throat_avg_bright) + mean(fem$throat_avg_bright)) 

ggplot(pred.fem.numEPO.throat.adj1) +
  geom_line(aes(x=throat_avg_bright_scaled, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=throat_avg_bright_scaled, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=throat_avg_bright_scaled, y=num.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Female throat average brightness") +
  ylab("Number of EPO")

ggplot(pred.fem.numEPO.throat.adj1) +
  geom_line(aes(x=throat.unscale, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=throat.unscale, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=throat_avg_bright, y=num.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Female throat average brightness") +
  ylab("Number of EPO")

# ggsave("output-files/fem num EPO model pred fem throat.png", h=4, w=4.5)


## Female breast color----------------------------------------------------------

# vector for range of breast
seq.fem.breast <- seq(range(fem$breast_avg_bright_scaled, na.rm=T)[1],
                      range(fem$breast_avg_bright_scaled, na.rm=T)[2], 0.001)
fem.breast.n <- length(seq.fem.breast)

# generate new data
pred.fem.numEPO.breast <-  
  data.frame(tail_scaled = rep(mean(fem$tail_scaled, narm=T), fem.breast.n),
             site.size = rep("medium", fem.breast.n),
       socM_tail_scaled = rep(mean(fem$socM_tail_scaled, na.rm=T), fem.breast.n),
       ci_1_julian_scaled = rep(mean(fem$ci_1_julian_scaled, na.rm=T), fem.breast.n),
       throat_avg_bright_scaled = rep(mean(fem$throat_avg_bright_scaled, na.rm=T), fem.breast.n),
       socM_t.avg.bright_scaled = rep(mean(fem$socM_t.avg.bright_scaled, na.rm=T), fem.breast.n),
       breast_avg_bright_scaled = seq.fem.breast,
       socM_r.avg.bright_scaled = rep(mean(fem$socM_r.avg.bright_scaled, na.rm=T), fem.breast.n))

# calculate fit and CIs
pred.fem.numEPO.breast.adj1 <- calculate.glm.ci(fem.numEPO.adj1,
                                                pred.fem.numEPO.breast)
# un-scale the predictor
pred.fem.numEPO.breast.adj1$breast.unscale <- (pred.fem.numEPO.breast.adj1$breast_avg_bright_scaled *
                                             sd(fem$breast_avg_bright, na.rm=T) + 
                                               mean(fem$breast_avg_bright, na.rm=T))

ggplot(pred.fem.numEPO.breast.adj1) +
  geom_line(aes(x=breast_avg_bright_scaled, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=breast_avg_bright_scaled, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=breast_avg_bright_scaled, y=num.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Female breast average brightness") +
  ylab("Number of EPO")

ggplot(pred.fem.numEPO.breast.adj1) +
  geom_line(aes(x=breast.unscale, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=breast.unscale, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=breast_avg_bright, y=num.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Female breast average brightness") +
  ylab("Number of EPO")


# ggsave("output-files/fem num EPO model pred fem breast.png", h=4, w=4.5)


# Visualize results for social mate --------------------------------------------

## Social male tail length

# vector for range of tail lengths
range(fem$socM_tail, na.rm=T) # 74, 103
seq.fem.Mtail <- seq(74, 103, 0.01)

seq.fem.Mtail <- seq(range(fem$socM_tail_scaled, na.rm=T)[1],
                      range(fem$socM_tail_scaled, na.rm=T)[2], 0.001)
fem.Mtail.n <- length(seq.fem.Mtail)

# generate new data
pred.fem.numEPO.Mtail <- 
  data.frame(tail_scaled = rep(mean(fem$tail_scaled, na.rm=T), fem.Mtail.n),
             site.size = rep("medium", fem.Mtail.n),
             socM_tail_scaled = seq.fem.Mtail,
             ci_1_julian_scaled = rep(mean(fem$ci_1_julian_scaled, na.rm=T), fem.Mtail.n),
             throat_avg_bright_scaled = rep(mean(fem$throat_avg_bright_scaled, na.rm=T), fem.Mtail.n),
             socM_t.avg.bright_scaled = rep(mean(fem$socM_t.avg.bright_scaled, na.rm=T), fem.Mtail.n),
             breast_avg_bright_scaled = rep(mean(fem$breast_avg_bright_scaled, na.rm=T), fem.Mtail.n),
             socM_r.avg.bright_scaled = rep(mean(fem$socM_r.avg.bright_scaled, na.rm=T), fem.Mtail.n))

# calculate fit and CIs
pred.fem.numEPO.Mtail.adj1 <- calculate.glm.ci(fem.numEPO.adj1, 
                                              pred.fem.numEPO.Mtail)
# un-scale the predictor
pred.fem.numEPO.Mtail.adj1$Mtail.unscale <- (pred.fem.numEPO.Mtail.adj1$socM_tail_scaled *
                                                 sd(fem$socM_tail, na.rm=T) + 
                                                 mean(fem$socM_tail, na.rm=T))

ggplot(pred.fem.numEPO.Mtail.adj1) +
  geom_line(aes(x=socM_tail_scaled, y=fit), size=1.5, linetype="dashed", color="darkgreen") +
  geom_ribbon(aes(x=socM_tail_scaled, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=socM_tail_scaled, y=num.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social male tail length") +
  ylab("Female number of EPO")

ggplot(pred.fem.numEPO.Mtail.adj1) +
  geom_line(aes(x=Mtail.unscale, y=fit), size=1.5, linetype="dashed", color="darkgreen") +
  geom_ribbon(aes(x=Mtail.unscale, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=socM_tail, y=num.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social male tail length") +
  ylab("Female number of EPO")

# ggsave("output-files/fem num EPO model pred fem Mtail.png", h=4, w=4.5)


## Social male throat-----------------------------------------------------------

# vector for range of throat color
seq.fem.Mthroat <- seq(range(fem$socM_t.avg.bright_scaled, na.rm=T)[1],
                       range(fem$socM_t.avg.bright_scaled, na.rm=T)[2], 0.001)
fem.Mthroat.n <- length(seq.fem.Mthroat)

# generate y-values
pred.fem.numEPO.Mthroat <-  
  data.frame(tail_scaled = rep(mean(fem$tail_scaled, narm=T), fem.Mthroat.n),
             site.size = rep("medium", fem.Mthroat.n),
             socM_tail_scaled = rep(mean(fem$socM_tail_scaled, na.rm=T), fem.Mthroat.n),
             ci_1_julian_scaled = rep(mean(fem$ci_1_julian_scaled, na.rm=T), fem.Mthroat.n),
             throat_avg_bright_scaled = rep(mean(fem$throat_avg_bright_scaled, na.rm=T), fem.Mthroat.n),
             socM_t.avg.bright_scaled = seq.fem.Mthroat,
             breast_avg_bright_scaled = rep(mean(fem$breast_avg_bright_scaled, na.rm=T), fem.Mthroat.n),
             socM_r.avg.bright_scaled = rep(mean(fem$socM_r.avg.bright_scaled, na.rm=T), fem.Mthroat.n))

# calculate fit and CIs
pred.fem.numEPO.Mthroat.adj1 <- calculate.glm.ci(fem.numEPO.adj1, 
                                                pred.fem.numEPO.Mthroat)
# un-scale the predictor
pred.fem.numEPO.Mthroat.adj1$Mthroat.unscale <- 
  (pred.fem.numEPO.Mthroat.adj1$socM_t.avg.bright_scaled *
                                               sd(fem$socM_t.avg.bright, na.rm=T) + 
                                               mean(fem$socM_t.avg.bright, na.rm=T))

ggplot(pred.fem.numEPO.Mthroat.adj1) +
  geom_line(aes(x=socM_t.avg.bright_scaled, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=socM_t.avg.bright_scaled, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=socM_t.avg.bright_scaled, y=num.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social male throat average brightness") +
  ylab("Female number of EPO")

ggplot(pred.fem.numEPO.Mthroat.adj1) +
  geom_line(aes(x=Mthroat.unscale, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=Mthroat.unscale, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=socM_t.avg.bright, y=num.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social male throat average brightness") +
  ylab("Female number of EPO")

# ggsave("output-files/fem num EPO model pred fem Mthroat.png", h=4, w=4.5)


## Social male breast----------------------------------------------------------- 

# vector for range of breast color
seq.fem.Mbreast <- seq(range(fem$socM_r.avg.bright_scaled, na.rm=T)[1],
                       range(fem$socM_r.avg.bright_scaled, na.rm=T)[2], 0.001)
fem.Mbreast.n <- length(seq.fem.Mbreast)

# generate y-values
pred.fem.numEPO.Mbreast <-  
  data.frame(tail_scaled = rep(mean(fem$tail_scaled, narm=T), fem.Mbreast.n),
             site.size = rep("medium", fem.Mbreast.n),
             socM_tail_scaled = rep(mean(fem$socM_tail_scaled, na.rm=T), fem.Mbreast.n),
             ci_1_julian_scaled = rep(mean(fem$ci_1_julian_scaled, na.rm=T), fem.Mbreast.n),
             throat_avg_bright_scaled = rep(mean(fem$throat_avg_bright_scaled, na.rm=T), fem.Mbreast.n),
             socM_t.avg.bright_scaled = rep(mean(fem$socM_t.avg.bright_scaled, na.rm=T), fem.Mbreast.n),
             breast_avg_bright_scaled = rep(mean(fem$breast_avg_bright_scaled, na.rm=T), fem.Mbreast.n),
             socM_r.avg.bright_scaled = seq.fem.Mbreast)

# calculate fit and CIs
pred.fem.numEPO.Mbreast.adj1 <- calculate.glm.ci(fem.numEPO.adj1, 
                                                 pred.fem.numEPO.Mbreast)
# un-scale the predictor
pred.fem.numEPO.Mbreast.adj1$Mbreast.unscale <- 
  (pred.fem.numEPO.Mbreast.adj1$socM_r.avg.bright_scaled *
     sd(fem$socM_r.avg.bright, na.rm=T) + 
     mean(fem$socM_r.avg.bright, na.rm=T))

ggplot(pred.fem.numEPO.Mbreast.adj1) +
  geom_line(aes(x=socM_r.avg.bright_scaled, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=socM_r.avg.bright_scaled, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=socM_r.avg.bright_scaled, y=num.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social male breast average brightness") +
  ylab("Female number of EPO")

ggplot(pred.fem.numEPO.Mbreast.adj1) +
  geom_line(aes(x=Mbreast.unscale, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=Mbreast.unscale, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=fem, aes(x=socM_r.avg.bright, y=num.epo), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social male breast average brightness") +
  ylab("Female number of EPO")

# ggsave("output-files/fem num EPO model pred fem Mbreast.png", h=4, w=4.5)

################################################################################

#-------------------------------------------------------------------------------
# Fit male model number EPO
#-------------------------------------------------------------------------------

# female throat was the only significant predictor in univariate models

male.numEPO.adj1 <- glm.nb(ep.chick ~ tail_scaled + socF_tail_scaled +
                             throat_avg_bright_scaled + socF_t.avg.bright_scaled + 
                             breast_avg_bright_scaled + socF_r.avg.bright_scaled +
                             ci_julian_scaled, data=male)

male.numEPO.adj2 <- glm.nb(ep.chick ~ site.size + tail_scaled + socF_tail_scaled +
                             throat_avg_bright_scaled + socF_t.avg.bright_scaled + 
                             breast_avg_bright_scaled + socF_r.avg.bright_scaled +
                             ci_julian_scaled, data=male, 
                           control=glm.control(maxit = 100))

# compare models with and without site
# support for dropping site
anova(male.numEPO.adj1, male.numEPO.adj2)


# get table of coefficients
table.male.numEPO.adj1 <- summary(male.numEPO.adj1)
table2.male.numEPO.adj1 <- as.data.frame(cbind(table.male.numEPO.adj1$coefficients, 
                              confint(male.numEPO.adj1)))
# back-transform
table2.male.numEPO.adj1$beta.BT <- exp(table2.male.numEPO.adj1$Estimate)
table2.male.numEPO.adj1$lowerCI.BT <- exp(table2.male.numEPO.adj1$`2.5 %`)
table2.male.numEPO.adj1$upperCI.BT <- exp(table2.male.numEPO.adj1$`97.5 %`)

write.csv(table2.male.numEPO.adj1, "output-files/male numEPO final mod table scaled.csv")

# get table of coefficients
table.male.numEPO.adj2 <- summary(male.numEPO.adj2)
table2.male.numEPO.adj2 <- as.data.frame(cbind(table.male.numEPO.adj2$coefficients, 
                                 confint(male.numEPO.adj2)))
# back-transform
table2.male.numEPO.adj2$beta.BT <- exp(table2.male.numEPO.adj2$Estimate)
table2.male.numEPO.adj2$lowerCI.BT <- exp(table2.male.numEPO.adj2$`2.5 %`)
table2.male.numEPO.adj2$upperCI.BT <- exp(table2.male.numEPO.adj2$`97.5 %`)

write.csv(table2.male.numEPO.adj2, "output-files/male numEPO final mod table site scaled.csv")


## Model without CI
male.numEPO.adj3 <- glm.nb(ep.chick ~ tail_scaled + socF_tail_scaled +
                             throat_avg_bright_scaled + socF_t.avg.bright_scaled + 
                             breast_avg_bright_scaled + socF_r.avg.bright_scaled,
                           data=male)
# get table of coefficients
table.male.numEPO.adj3 <- summary(male.numEPO.adj3)
table2.male.numEPO.adj3 <- as.data.frame(cbind(table.male.numEPO.adj3$coefficients, 
                                 confint(male.numEPO.adj3)))
# back-transform
table2.male.numEPO.adj3$beta.BT <- exp(table2.male.numEPO.adj3$Estimate)
table2.male.numEPO.adj3$lowerCI.BT <- exp(table2.male.numEPO.adj3$`2.5 %`)
table2.male.numEPO.adj3$upperCI.BT <- exp(table2.male.numEPO.adj3$`97.5 %`)


write.csv(table2.male.numEPO.adj3, "output-files/male numEPO final table no ci scaled.csv")

###-----------------------------------------------------------------------------
### visualize results for male models
###-----------------------------------------------------------------------------

## social female throat color---------------------------------------------------

# range
seq.male.Fthroat <- seq(range(male$socF_t.avg.bright_scaled)[1],
                        range(male$socF_t.avg.bright_scaled)[2], 0.001)
male.Fthroat.n <- length(seq.male.Fthroat)

# generate new data
pred.male.numEPO.Fthroat <- 
  data.frame(socF_t.avg.bright_scaled = seq.male.Fthroat,
             site.size = rep("medium",male.Fthroat.n ),
       throat_avg_bright_scaled = rep(mean(male$throat_avg_bright_scaled, na.rm=T), male.Fthroat.n),
       ci_julian_scaled = rep(mean(male$ci_julian_scaled , na.rm=T), male.Fthroat.n),
       tail_scaled = rep(mean(male$tail_scaled, na.rm=T), male.Fthroat.n),
       socF_tail_scaled = rep(mean(male$socF_tail_scaled, na.rm=T), male.Fthroat.n),
       breast_avg_bright_scaled = rep(mean(male$breast_avg_bright_scaled, na.rm=T), male.Fthroat.n),
       socF_r.avg.bright_scaled = rep(mean(male$socF_r.avg.bright_scaled, na.rm=T), male.Fthroat.n))

# calculate fit and CIs
pred.male.numEPO.Fthroat.adj1 <- calculate.glm.ci(male.numEPO.adj1, 
                                                 pred.male.numEPO.Fthroat)
# un-scale the predictor
pred.male.numEPO.Fthroat.adj1$Fthroat.unscale <- 
  (pred.male.numEPO.Fthroat.adj1$socF_t.avg.bright_scaled *
     sd(male$socF_t.avg.bright) + mean(male$socF_t.avg.bright))

ggplot(pred.male.numEPO.Fthroat.adj1) +
  geom_line(aes(x=socF_t.avg.bright_scaled, y=fit), size=1.5, linetype="solid",
            color="darkgreen") +
  geom_ribbon(aes(x=socF_t.avg.bright_scaled, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=socF_t.avg.bright_scaled, y=ep.chick), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social female throat average brightness") +
  ylab("Male number of EPO")

ggplot(pred.male.numEPO.Fthroat.adj1) +
  geom_line(aes(x=Fthroat.unscale, y=fit), size=1.5, linetype="solid",
            color="darkgreen") +
  geom_ribbon(aes(x=Fthroat.unscale, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=socF_t.avg.bright, y=ep.chick), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social female throat average brightness") +
  ylab("Male number of EPO")


# ggsave("output-files/male num EPO model pred socF throat.png", h=4, w=4.5)


## social female tail-----------------------------------------------------------

# range
seq.male.Ftail <- seq(range(male$socF_tail_scaled)[1],
                      range(male$socF_tail_scaled)[2], 0.001)
male.Ftail.n <- length(seq.male.Ftail)

# generate new data
pred.male.numEPO.Ftail <- 
  data.frame(socF_t.avg.bright_scaled = rep(mean(male$socF_t.avg.bright_scaled, na.rm=T), male.Ftail.n),
             site.size = rep("medium",male.Ftail.n ),
             throat_avg_bright_scaled = rep(mean(male$throat_avg_bright_scaled, na.rm=T), male.Ftail.n),
             ci_julian_scaled = rep(mean(male$ci_julian_scaled , na.rm=T), male.Ftail.n),
             tail_scaled = rep(mean(male$tail_scaled, na.rm=T), male.Ftail.n),
             socF_tail_scaled = seq.male.Ftail,
             breast_avg_bright_scaled = rep(mean(male$breast_avg_bright_scaled, na.rm=T), male.Ftail.n),
             socF_r.avg.bright_scaled = rep(mean(male$socF_r.avg.bright_scaled, na.rm=T), male.Ftail.n))

# calculate fit and CIs
pred.male.numEPO.Ftail.adj1 <- calculate.glm.ci(male.numEPO.adj1, 
                                                  pred.male.numEPO.Ftail)
# un-scale the predictor
pred.male.numEPO.Ftail.adj1$Ftail.unscale <- 
  (pred.male.numEPO.Ftail.adj1$socF_tail_scaled *
     sd(male$socF_tail) + mean(male$socF_tail))

ggplot(pred.male.numEPO.Ftail.adj1) +
  geom_line(aes(x=socF_tail_scaled, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=socF_tail_scaled, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=socF_tail_scaled, y=ep.chick), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social female tail streamer length") +
  ylab("Male number of EPO")

ggplot(pred.male.numEPO.Ftail.adj1) +
  geom_line(aes(x=Ftail.unscale, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=Ftail.unscale, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=socF_tail, y=ep.chick), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social female tail streamer length") +
  ylab("Male number of EPO")


# ggsave("output-files/male num EPO model pred socF tail.png", h=4, w=4.5)


## Social female breast---------------------------------------------------------

# range
seq.male.Fbreast <- seq(range(male$socF_r.avg.bright_scaled, na.rm=T)[1],
                        range(male$socF_r.avg.bright_scaled, na.rm=T)[2], 0.001)
male.Fbreast.n <- length(seq.male.Fbreast)

# generate new data
pred.male.numEPO.Fbreast <- 
  data.frame(socF_t.avg.bright_scaled = rep(mean(male$socF_t.avg.bright_scaled, na.rm=T), 
                                            male.Fbreast.n),
             site.size = rep("medium",male.Fbreast.n ),
             throat_avg_bright_scaled = rep(mean(male$throat_avg_bright_scaled, na.rm=T), 
                                            male.Fbreast.n),
             ci_julian_scaled = rep(mean(male$ci_julian_scaled , na.rm=T), male.Fbreast.n),
             tail_scaled = rep(mean(male$tail_scaled, na.rm=T), male.Fbreast.n),
             socF_tail_scaled = rep(mean(male$socF_tail_scaled, na.rm=T), male.Fbreast.n),
             breast_avg_bright_scaled = rep(mean(male$breast_avg_bright_scaled, na.rm=T), 
                                            male.Fbreast.n),
             socF_r.avg.bright_scaled = seq.male.Fbreast)

# calculate fit and CIs
pred.male.numEPO.Fbreast.adj1 <- calculate.glm.ci(male.numEPO.adj1, 
                                                  pred.male.numEPO.Fbreast)
# un-scale the predictor
pred.male.numEPO.Fbreast.adj1$Fbreast.unscale <- 
  (pred.male.numEPO.Fbreast.adj1$socF_r.avg.bright_scaled *
     sd(male$socF_r.avg.bright, na.rm=T) + mean(male$socF_r.avg.bright, na.rm=T))

ggplot(pred.male.numEPO.Fbreast.adj1) +
  geom_line(aes(x=socF_r.avg.bright_scaled, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=socF_r.avg.bright_scaled, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=socF_r.avg.bright_scaled, y=ep.chick), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social female breast average brightness") +
  ylab("Male number of EPO")

ggplot(pred.male.numEPO.Fbreast.adj1) +
  geom_line(aes(x=Fbreast.unscale, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=Fbreast.unscale, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=socF_r.avg.bright, y=ep.chick), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Social female breast average brightness") +
  ylab("Male number of EPO")


# ggsave("output-files/male num EPO model pred socF breast.png", h=4, w=4.5)


## male tail--------------------------------------------------------------------

# range of lengths
seq.male.tail <- seq(range(male$tail_scaled, na.rm=T)[1],
                     range(male$tail_scaled, na.rm=T)[2], 0.001)
male.tail.n <- length(seq.male.tail)

# generate new data
pred.male.numEPO.tail <- data.frame(tail_scaled = seq.male.tail,
                                    site.size = rep("medium", male.tail.n),
       socF_tail_scaled = rep(mean(male$socF_tail_scaled, na.rm=T), male.tail.n),
       ci_julian_scaled = rep(mean(male$ci_julian_scaled, na.rm=T), male.tail.n),
       throat_avg_bright_scaled = rep(mean(male$throat_avg_bright_scaled, na.rm=T), male.tail.n),
       socF_t.avg.bright_scaled = rep(mean(male$socF_t.avg.bright_scaled, na.rm=T), male.tail.n),
       breast_avg_bright_scaled = rep(mean(male$breast_avg_bright_scaled, na.rm=T), male.tail.n),
       socF_r.avg.bright_scaled = rep(mean(male$socF_r.avg.bright_scaled, na.rm=T), male.tail.n))

# calculate CI and fit
pred.male.numEPO.tail.adj1 <- calculate.glm.ci(male.numEPO.adj1,
                                               pred.male.numEPO.tail)
# un-scale the predictor
pred.male.numEPO.tail.adj1$tail.unscale <- 
  (pred.male.numEPO.tail.adj1$tail_scaled *
     sd(male$tail, na.rm=T) + mean(male$tail, na.rm=T))

ggplot(pred.male.numEPO.tail.adj1) +
  geom_line(aes(x=tail_scaled, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=tail_scaled, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=tail_scaled, y=ep.chick), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Male tail length") +
  ylab("Number of EPO")

ggplot(pred.male.numEPO.tail.adj1) +
  geom_line(aes(x=tail.unscale, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=tail.unscale, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=tail, y=ep.chick), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Male tail length") +
  ylab("Number of EPO")


# ggsave("output-files/male num EPO model pred male tail.png", h=4, w=4.5)


## Male throat------------------------------------------------------------------ 

# range
seq.male.throat <- seq(range(male$throat_avg_bright_scaled, na.rm=T)[1],
                       range(male$throat_avg_bright_scaled, na.rm=T)[2], 0.001)
male.throat.n <- length(seq.male.throat)

# generate y-values
pred.male.numEPO.throat <- data.frame(tail_scaled = rep(mean(male$tail_scaled, na.rm=T), male.throat.n),
                                      site.size=rep("medium", male.throat.n),
       socF_tail_scaled = rep(mean(male$socF_tail_scaled, na.rm=T), male.throat.n),
       ci_julian_scaled = rep(mean(male$ci_julian_scaled, na.rm=T), male.throat.n),
       throat_avg_bright_scaled = seq.male.throat,
       socF_t.avg.bright_scaled = rep(mean(male$socF_t.avg.bright_scaled, na.rm=T), male.throat.n),
       breast_avg_bright_scaled = rep(mean(male$breast_avg_bright_scaled, na.rm=T), male.throat.n),
       socF_r.avg.bright_scaled = rep(mean(male$socF_r.avg.bright_scaled, na.rm=T), male.throat.n))

# calculate fit and CIs
pred.male.numEPO.throat.adj1 <- calculate.glm.ci(male.numEPO.adj1,
                                                 pred.male.numEPO.throat)
# un-scale the predictor
pred.male.numEPO.throat.adj1$throat.unscale <- 
  (pred.male.numEPO.throat.adj1$throat_avg_bright_scaled *
     sd(male$throat_avg_bright, na.rm=T) + mean(male$throat_avg_bright, na.rm=T))

ggplot(pred.male.numEPO.throat.adj1) +
  geom_line(aes(x=throat_avg_bright_scaled, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=throat_avg_bright_scaled, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=throat_avg_bright_scaled, y=ep.chick), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Male throat average brightness") +
  ylab("Male number of EPO")

ggplot(pred.male.numEPO.throat.adj1) +
  geom_line(aes(x=throat.unscale, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=throat.unscale, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=throat_avg_bright, y=ep.chick), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Male throat average brightness") +
  ylab("Male number of EPO")

# ggsave("output-files/male num EPO model pred male throat.png", h=4, w=4.5)


## Male breast------------------------------------------------------------------

# range
seq.male.breast <- seq(range(male$breast_avg_bright_scaled, na.rm=T)[1],
                       range(male$breast_avg_bright_scaled, na.rm=T)[2], 0.001)
male.breast.n <- length(seq.male.breast)

# generate y-values
pred.male.numEPO.breast <- data.frame(tail_scaled = rep(mean(male$tail_scaled, na.rm=T), male.breast.n),
                                      site.size = rep("medium", male.breast.n),
       socF_tail_scaled = rep(mean(male$socF_tail_scaled, na.rm=T), male.breast.n),
       ci_julian_scaled = rep(mean(male$ci_julian_scaled, na.rm=T), male.breast.n),
       throat_avg_bright_scaled = rep(mean(male$throat_avg_bright_scaled, na.rm=T), male.breast.n),
       socF_t.avg.bright_scaled = rep(mean(male$socF_t.avg.bright_scaled, na.rm=T), male.breast.n),
       breast_avg_bright_scaled = seq.male.breast,
       socF_r.avg.bright_scaled = rep(mean(male$socF_r.avg.bright_scaled, na.rm=T), male.breast.n))

# calculate fit and CIs
pred.male.numEPO.breast.adj1 <- calculate.glm.ci(male.numEPO.adj1, 
                                                 pred.male.numEPO.breast)
# un-scale the predictor
pred.male.numEPO.breast.adj1$breast.unscale <- 
  (pred.male.numEPO.breast.adj1$breast_avg_bright_scaled *
     sd(male$breast_avg_bright, na.rm=T) + mean(male$breast_avg_bright, na.rm=T))

ggplot(pred.male.numEPO.breast.adj1) +
  geom_line(aes(x=breast_avg_bright_scaled, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=breast_avg_bright_scaled, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=breast_avg_bright_scaled, y=ep.chick), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Male breast average brightness") +
  ylab("Male number of EPO")

ggplot(pred.male.numEPO.breast.adj1) +
  geom_line(aes(x=breast.unscale, y=fit), size=1.5, linetype="dashed",
            color="darkgreen") +
  geom_ribbon(aes(x=breast.unscale, ymax=upper, ymin=lower), 
              fill="darkgreen", alpha=0.3) +
  geom_point(data=male, aes(x=breast_avg_bright, y=ep.chick), 
             alpha=0.5, size=3, color="darkgreen") +
  xlab("Male breast average brightness") +
  ylab("Male number of EPO")

# ggsave("output-files/male num EPO model pred male breast.png", h=4, w=4.5)


#-------------------------------------------------------------------------------
# Visualize female and male num EPO together
#-------------------------------------------------------------------------------

# set colors to create a legend
cols <- c("Females"="#CC00FF", "Males"="#00B0F0")


## Tail ------------------------------------------------------------------------

numEPO.focal.tail <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.numEPO.tail.adj1, 
            aes(x=tail.unscale, y=fit, color="Females"), 
            size=1.5, linetype="solid") +
  geom_ribbon(data=pred.fem.numEPO.tail.adj1, 
              aes(x=tail.unscale, ymin=lower, ymax=upper, fill="Females"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=tail, y=num.epo, color="Females"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.numEPO.tail.adj1, 
            aes(x=tail.unscale, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.numEPO.tail.adj1, 
              aes(x=tail.unscale, ymax=upper, ymin=lower, fill="Males"),
              alpha=0.3) +
  geom_point(data=male, aes(x=tail, y=ep.chick, color="Males"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL) +
  theme_light()

ggsave("output-files/combined plot female and male tail on num EPO.png", 
       h=4, w=7, scale=0.5)



# Throat -----------------------------------------------------------------------

numEPO.focal.throat <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.numEPO.throat.adj1, 
            aes(x=throat.unscale, y=fit, color="Females"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.numEPO.throat.adj1, 
              aes(x=throat.unscale, ymin=lower, ymax=upper, fill="Females"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=throat_avg_bright, y=num.epo, color="Females"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.numEPO.throat.adj1, 
            aes(x=throat.unscale, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.numEPO.throat.adj1, 
              aes(x=throat.unscale, ymax=upper, ymin=lower, fill="Males"),
              alpha=0.3) +
  geom_point(data=male, aes(x=throat_avg_bright, y=ep.chick, color="Males"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL) +
  theme_light()


ggsave("output-files/combined plot female and male throat on num EPO.png", 
       h=4, w=7, scale=0.5)

# Breast -----------------------------------------------------------------------

numEPO.focal.breast <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.numEPO.breast.adj1, 
            aes(x=breast.unscale, y=fit, color="Females"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.numEPO.breast.adj1, 
              aes(x=breast.unscale, ymin=lower, ymax=upper, fill="Females"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=breast_avg_bright, y=num.epo, color="Females"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.numEPO.breast.adj1, 
            aes(x=breast.unscale, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.numEPO.breast.adj1, 
              aes(x=breast.unscale, ymax=upper, ymin=lower, fill="Males"),
              alpha=0.3) +
  geom_point(data=male, aes(x=breast_avg_bright, y=ep.chick, color="Males"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL) +
  theme_light()


ggsave("output-files/combined plot female and male breast on num EPO.png", 
       h=4, w=7, scale=0.5)



## social mate tail-------------------------------------------------------------

numEPO.social.tail <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.numEPO.Mtail.adj1, 
            aes(x=Mtail.unscale, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.numEPO.Mtail.adj1, 
              aes(x=Mtail.unscale, ymin=lower, ymax=upper, fill="Males"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=socM_tail, y=num.epo, color="Males"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.numEPO.Ftail.adj1, 
            aes(x=Ftail.unscale, y=fit, color="Females"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.numEPO.Ftail.adj1, 
              aes(x=Ftail.unscale, ymax=upper, ymin=lower, fill="Females"),
              alpha=0.3) +
  geom_point(data=male, aes(x=socF_tail, y=ep.chick, color="Females"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL)+
  theme_light()

ggsave("output-files/combined plot female and male mate tail on numEPO.png", 
       h=4, w=7, scale=0.5)


## social mate throat-----------------------------------------------------------

numEPO.social.throat <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.numEPO.Mthroat.adj1, 
            aes(x=Mthroat.unscale, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.numEPO.Mthroat.adj1, 
              aes(x=Mthroat.unscale, ymin=lower, ymax=upper, fill="Males"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=socM_t.avg.bright, y=num.epo, color="Males"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.numEPO.Fthroat.adj1, 
            aes(x=Fthroat.unscale, y=fit, color="Females"), 
            size=1.5, linetype="solid") +
  geom_ribbon(data=pred.male.numEPO.Fthroat.adj1, 
              aes(x=Fthroat.unscale, ymax=upper, ymin=lower, fill="Females"),
              alpha=0.3) +
  geom_point(data=male, aes(x=socF_t.avg.bright, y=ep.chick, color="Females"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL)+
  theme_light()

ggsave("output-files/combined plot female and male mate throat on numEPO.png", 
       h=4, w=7, scale=0.5)


## social mate breast-----------------------------------------------------------

numEPO.social.breast <- ggplot() +
  # female predictions and points
  geom_line(data=pred.fem.numEPO.Mbreast.adj1, 
            aes(x=Mbreast.unscale, y=fit, color="Males"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.fem.numEPO.Mbreast.adj1, 
              aes(x=Mbreast.unscale, ymin=lower, ymax=upper, fill="Males"),
              alpha=0.3) +
  geom_point(data=fem, aes(x=socM_r.avg.bright, y=num.epo, color="Males"),
             alpha=0.5) +
  # male predictions and points             
  geom_line(data=pred.male.numEPO.Fbreast.adj1, 
            aes(x=Fbreast.unscale, y=fit, color="Females"), 
            size=1.5, linetype="dashed") +
  geom_ribbon(data=pred.male.numEPO.Fbreast.adj1, 
              aes(x=Fbreast.unscale, ymax=upper, ymin=lower, fill="Females"),
              alpha=0.3) +
  geom_point(data=male, aes(x=socF_r.avg.bright, y=ep.chick, color="Females"),
             alpha=0.5) +
  scale_color_manual(name="Legend", values=cols) +
  scale_fill_manual(name="Legend", values=cols) +
  xlab(NULL) +
  ylab(NULL)+
  theme_light()

ggsave("output-files/combined plot female and male mate breast on numEPO.png", 
       h=4, w=7, scale=0.5)

## Save plots to load in a different script

save(numEPO.focal.breast, numEPO.focal.tail, numEPO.focal.throat, 
     numEPO.social.breast, numEPO.social.tail, numEPO.social.throat,
     file="output-files/numEPO combined plots_2025-03-13.Rdata")



