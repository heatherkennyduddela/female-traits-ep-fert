
################################################################################
# Script for making combined figure panels of EP mating results
# Heather Kenny-Duddela
# Feb 11, 2025
################################################################################

# load data

load("output-files/numEPO combined plots_2025-10-02_bluetan.Rdata") # numEPO plots
load("output-files/numMates combined plots_2025-03-13_bluetan.Rdata") # number Mates plots

# libraries
library(ggplot2) # for plotting
library(ggpubr) # for arranging multiple plots in one panel

#-------------------------------------------------------------------------------

### Panel for focal bird trait effects

# remove legends and add axis labels for some


# Num EPO
numEPO.focal.tail2 <- numEPO.focal.tail + guides(color="none", fill="none") +
  ylab("Number of EPO") + ggtitle(" ")

numEPO.focal.throat2 <- numEPO.focal.throat + guides(color="none", fill="none") +
  ggtitle(" ")

numEPO.focal.breast2 <- numEPO.focal.breast + guides(color="none", fill="none") +
  ggtitle(" ")



# Num mates
focal.tail2 <- focal.tail + guides(color="none", fill="none") +
  ylab("Number of genetic mates") +
  xlab("Tail streamer length (mm)")

focal.throat2 <- focal.throat + guides(color="none", fill="none") +
  xlab("Throat average brightness")

focal.breast2 <- focal.breast + guides(color="none", fill="none") +
  xlab("Breast average brightness")

# combine panels
ggarrange(numEPO.focal.tail2, numEPO.focal.throat2, numEPO.focal.breast2,
          focal.tail2, focal.throat2, focal.breast2,
          labels=c("A","B","C","D","E","F"),
          nrow=2, ncol=3, align="hv",
          label.x=0.15, label.y=1)

ggsave("output-files/focal bird trait effects all panels_2025-10-02_bluetan.png",
       h=5, w=7)


#-------------------------------------------------------------------------------
# Panel for social mate trait effects

# remove legends and add axis labels for some


# Num EPO
numEPO.social.tail2 <- numEPO.social.tail + guides(color="none", fill="none") +
  ylab("Number of EPO") + ggtitle("    ")

numEPO.social.throat2 <- numEPO.social.throat + guides(color="none", fill="none") +
  ggtitle("    ")

numEPO.social.breast2 <- numEPO.social.breast + guides(color="none", fill="none") +
  ggtitle("    ")


# Num mates
social.tail2 <- social.tail + guides(color="none", fill="none") +
  ylab("Number of genetic mates") +
  xlab("Tail streamer length (mm)")

social.throat2 <- social.throat + guides(color="none", fill="none") +
  xlab("Throat average brightness")

social.breast2 <- social.breast + guides(color="none", fill="none") +
  xlab("Breast average brightness")

# combine panels
ggarrange(numEPO.social.tail2, numEPO.social.throat2, numEPO.social.breast2,
          social.tail2, social.throat2, social.breast2,
          labels=c("A","B","C","D","E","F"),
          nrow=2, ncol=3, align="hv",
          label.x=0.15, label.y=1)

ggsave("output-files/social mate trait effects all panels_2025-10-02_bluetan.png",
       h=5, w=7)






