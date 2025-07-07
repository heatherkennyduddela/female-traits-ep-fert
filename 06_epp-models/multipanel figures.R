
################################################################################
# Script for making combined figure panels of EP mating results
# Heather Kenny-Duddela
# Feb 11, 2025
################################################################################

# load data

load("output-files/binary EP combined plots_2025-03-13_bluetan.Rdata") # binary EP plots
load("output-files/numEPO combined plots_2025-03-13_bluetan.Rdata") # numEPO plots
load("output-files/numMates combined plots_2025-03-13_bluetan.Rdata") # number Mates plots
load("output-files/prop EP combined plots_2025-03-13_bluetan.Rdata") # proportion EPO plots

# libraries
library(ggplot2) # for plotting
library(ggpubr) # for arranging multiple plots in one panel

#-------------------------------------------------------------------------------

### Panel for focal bird trait effects

# remove legends and add axis labels for some

# binary EPP
binary.focal.tail2 <- binary.focal.tail + guides(color="none", fill="none") +
  ylab("Binary EP status") + ggtitle("    ")

binary.focal.throat2 <- binary.focal.throat + 
  guides(color="none", fill="none") + ggtitle("    ")

binary.focal.breast2 <- binary.focal.breast + 
  guides(color="none", fill="none") + ggtitle("    ")

# Num EPO
numEPO.focal.tail2 <- numEPO.focal.tail + guides(color="none", fill="none") +
  ylab("Number of EPO")

numEPO.focal.throat2 <- numEPO.focal.throat + guides(color="none", fill="none")

numEPO.focal.breast2 <- numEPO.focal.breast + guides(color="none", fill="none")

# Prop EPO
prop.focal.tail2 <- prop.focal.tail + guides(color="none", fill="none") +
  ylab("Proportion EPO")

prop.focal.throat2 <- prop.focal.throat + guides(color="none", fill="none")

prop.focal.breast2 <- prop.focal.breast + guides(color="none", fill="none")

# Num mates
focal.tail2 <- focal.tail + guides(color="none", fill="none") +
  ylab("Number of genetic mates") +
  xlab("Tail streamer length (mm)")

focal.throat2 <- focal.throat + guides(color="none", fill="none") +
  xlab("Throat average brightness")

focal.breast2 <- focal.breast + guides(color="none", fill="none") +
  xlab("Breast average brightness")

# combine panels
ggarrange(binary.focal.tail2, binary.focal.throat2, binary.focal.breast2,
          numEPO.focal.tail2, numEPO.focal.throat2, numEPO.focal.breast2,
          prop.focal.tail2, prop.focal.throat2, prop.focal.breast2,
          focal.tail2, focal.throat2, focal.breast2,
          labels=c("A","B","C","D","E","F","G","H","I","J","K","L"),
          nrow=4, ncol=3, align="hv",
          label.x=0.25, label.y=1)

ggsave("output-files/focal bird trait effects all panels_2025-03-13_bluetan.png",
       h=8, w=7)


#-------------------------------------------------------------------------------
# Panel for social mate trait effects

# remove legends and add axis labels for some

# binary EPP
binary.social.tail2 <- binary.social.tail + guides(color="none", fill="none") +
  ylab("Binary EP status") + ggtitle("    ")

binary.social.throat2 <- binary.social.throat + 
  guides(color="none", fill="none") + ggtitle("    ")

binary.social.breast2 <- binary.social.breast + 
  guides(color="none", fill="none") + ggtitle("    ")

# Num EPO
numEPO.social.tail2 <- numEPO.social.tail + guides(color="none", fill="none") +
  ylab("Number of EPO")

numEPO.social.throat2 <- numEPO.social.throat + guides(color="none", fill="none")

numEPO.social.breast2 <- numEPO.social.breast + guides(color="none", fill="none")

# Prop EPO
prop.social.tail2 <- prop.social.tail + guides(color="none", fill="none") +
  ylab("Proportion EPO")

prop.social.throat2 <- prop.social.throat + guides(color="none", fill="none")

prop.social.breast2 <- prop.social.breast + guides(color="none", fill="none")

# Num mates
social.tail2 <- social.tail + guides(color="none", fill="none") +
  ylab("Number of genetic mates") +
  xlab("Tail streamer length (mm)")

social.throat2 <- social.throat + guides(color="none", fill="none") +
  xlab("Throat average brightness")

social.breast2 <- social.breast + guides(color="none", fill="none") +
  xlab("Breast average brightness")

# combine panels
ggarrange(binary.social.tail2, binary.social.throat2, binary.social.breast2,
          numEPO.social.tail2, numEPO.social.throat2, numEPO.social.breast2,
          prop.social.tail2, prop.social.throat2, prop.social.breast2,
          social.tail2, social.throat2, social.breast2,
          labels=c("A","B","C","D","E","F","G","H","I","J","K","L"),
          nrow=4, ncol=3, align="hv",
          label.x=0.25, label.y=1)

ggsave("output-files/social mate trait effects all panels_2025-03-13_bluetan.png",
       h=8, w=7)






