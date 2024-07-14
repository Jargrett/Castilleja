
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU")

library(tidyverse)
library(lme4)
library(ggplot2)
library(ggpubr)
library(car)
library(vegan)

#Load in 2023 + 2024 datasets (For diversity analysis we will use species counts)
CALI.24 <- read.csv("Cali 2024 Count.csv")
CALI.23 <- read.csv("Cali 2023 Count.csv")

#Adding a Column for year
CALI.23 <- CALI.23 %>%
  mutate(Year = 2023)
  
CALI.24 <- CALI.24 %>%
  mutate(Year = 2024)

#sperating the species matrix from the environmental data
CALI.23.env <- subset(CALI.23, select=c(1,2,4:6,50)) #gathering enviornmental data
CALI.24.env <- subset(CALI.24, select=c(1,2,4:6,70))

CALI.23.species <- CALI.23[ -c(1:6,50)] #Isolating the species matrix
CALI.24.species <- CALI.24[ -c(1:6,70)]
  
#--------------------Diversity Analysis--------------------#
# calculating diversity metrics for our linariifolia sites

# Calculating Shannon diversity,richness, and evenness for 2023 plots
div.23 <- diversity(CALI.23.species, index = "shannon")
rich.23 <- specnumber(CALI.23.species)
even.23 <- diversity(CALI.23.species, index = "shannon") / log(specnumber(CALI.23.species)) 

# Calculating Shannon diversity,richness, and eveness for 2024 plots
div.24 <- diversity(CALI.24.species, index = "shannon")
rich.24 <- specnumber(CALI.24.species)
even.24 <- diversity(CALI.24.species, index = "shannon") / log(specnumber(CALI.24.species)) 

#combined data set with environmental and calculated values
CALI.23.div <- cbind(CALI.23.env,div.23,rich.23,even.23)
CALI.24.div <- cbind(CALI.24.env,div.24,rich.24,even.24)

CALI.23.div <- CALI.23.div %>% #renaming columns to prepare for merging
  rename("div" = "div.23",
         "even" = "even.23",
         "rich" = "rich.23")

CALI.24.div <- CALI.24.div %>% 
  rename("div" = "div.24",
         "even" = "even.24",
         "rich" = "rich.24")

#merging dataframes
cali <- rbind(CALI.23.div,CALI.24.div)
cali <- cali %>%
  relocate(Year)
cali <- cali %>% 
  rename("Castilleja" = "Treatment")

#creating a new csv for later use and sharing
write.csv(cali, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined linariifolia Diversity.csv", row.names=FALSE)

#--------Models------#
cali.div <- lmer(div ~ Castilleja*Site + (1|Pair) + (1|Year), data = cali)
summary(cali.div)
Anova(cali.div)

cali.rich <- lmer(rich ~ Castilleja*Site + (1|Pair) + (1|Year), data = cali)
summary(cali.rich)
Anova(cali.rich)

cali.even <- lmer(even ~ Castilleja*Site + (1|Pair) + (1|Year), data = cali)
summary(cali.even)
Anova(cali.even)

#-------plotting------#
cali.diversity <- ggplot(cali, aes(x = Castilleja, y = div)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~Site) +
  scale_color_manual(values=c( "coral", "burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Shannon Diversity") +
  ylim(0,3)

cali.diversity


cali.richness <- ggplot(cali, aes(x = Castilleja, y = rich)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~Site) +
  scale_color_manual(values=c( "coral", "burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Species Richness") +
  ylim(0,20)

cali.richness

cali.evenness <- ggplot(cali, aes(x = Castilleja, y = even)) +
  stat_summary(aes(group = Pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (Castilleja), size = 1, alpha = 2), show.legend = FALSE) + 
  #stat_summary(fun.y=mean, position = position_dodge(1), size = 2, col = "grey34") +
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~Site) +
  scale_color_manual(values=c( "coral", "burlywood4")) +
  labs(x = "Castilleja linariifolia", y = "Species Evenness") +
  ylim(.4,1)

cali.evenness

linariifolia.plot <- ggarrange(cali.diversity, cali.richness, cali.evenness,
                         labels = c("A", "B","C"), 
                         nrow = 1, common.legend = TRUE, legend = "bottom")
linariifolia.plot



