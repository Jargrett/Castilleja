# Almont Chromosa analysis
# Initiated: 6/26/24
# Completed: TBD

#####################
#                   # - Impact of Castilleja presence on:
#      Question     # - richness and evenness (Shannon Diversity Index)
#                   # - 
##################### - Also would be good too look at: Individual species cover

# Set working directory (Workspace)
# This can be done manually or through session -> set working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Anders Data 2024")

# Load-in packages
# I will explain these if needed
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(tidyverse) # for data working
library(vegan) # for diversity
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car)
library(lme4)

# importing dataframe

al <- read.csv("CACH Estimated Cover.csv")
al.cover <- al[ -c(1:5,7,56:60)] # removing unnecessary columns 
al.env <- subset(al, select=c(2,3,5:7)) # holding our environmental stuff

#--------------------Diversity Analysis--------------------#
# calculating diversity metrics for our linariifolia sites

# Calculating the species richness for plots
rich.cov <- specnumber(al.cover)

# Calculating Shannon diversity for plots
div.cov <- diversity(al.cover, index = "shannon")


# Calculating species evenness for plots
even.cov <- diversity(al.cover, index = "shannon") / log(specnumber(al.cover))  

#combined data set with environmental and calculated values
al.div <- cbind(al.env,rich.cov,div.cov,even.cov)

#removing NA rows
al.cov.div <- subset(al.div,!is.na(Site))

#renaming columns for stats later
al.cov.div <- al.cov.div %>% rename(Castilleja = Castilleja.)

#export this dataset for combined analysis
write.csv(al.cov.div, "C:\\Users\\jargr\\Dropbox\\PC\\Desktop\\Data Analysis\\RMBL\\Anders Data 2024\\Almont 2024 Diversity.csv")


#--------Models------#
cach.div <- lmer(div.cov ~ Castilleja*Site + (1|Paired.Plot), data = al.cov.div)
summary(cach.div)
Anova(cach.div)

cach.rich <- lmer(rich.cov ~ Castilleja*Site + (1|Paired.Plot), data = al.cov.div)
summary(cach.rich)
Anova(cach.rich)

cach.even <- lmer(even.cov ~ Castilleja*Site + (1|Paired.Plot), data = al.cov.div)
summary(cach.even)
Anova(cach.even)

#----Diversity plotting----#

almont.div <- ggplot(al.cov.div, aes(x = Castilleja, y = div.cov)) +
  stat_summary(aes(group = Paired.Plot), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  scale_color_manual(values=c("burlywood4", "coral3")) +
  labs(x = "Castilleja chromosa", y = "Shannon Diversity") +
  ylim(1,2.5)

almont.div

almont.rich <- ggplot(al.cov.div, aes(x = Castilleja, y = rich.cov)) +
  stat_summary(aes(group = Paired.Plot), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  scale_color_manual(values=c("burlywood4", "coral3")) +
  labs(x = "Castilleja chromosa", y = "Species Richness") +
  ylim(8,20)

almont.rich

almont.even <- ggplot(al.cov.div, aes(x = Castilleja, y = even.cov)) +
  stat_summary(aes(group = Paired.Plot), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  scale_color_manual(values=c("burlywood4", "coral3")) +
  labs(x = "Castilleja chromosa", y = "Species Evenness") +
  ylim(.4,1)

almont.even

almont.plot <- ggarrange(almont.div, almont.rich, almont.even,
                         labels = c("A", "B","C"), 
                         nrow = 1, common.legend = TRUE, legend = "bottom")
almont.plot


#----------Plant Species Cover analysis-----------#
