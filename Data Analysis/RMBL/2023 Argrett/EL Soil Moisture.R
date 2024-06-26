# EL Soil Moisture Calulations
# Initiated: 8/23/2023
# Complted:

setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2023 Argrett")

# Load-in packages
library(tidyverse) # for data working
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car) # Regression Analysis

soil.moisture <- read.csv("Emerald Lake Soil Moisture Data.csv")

#checking data structure
str(soil.moisture)
soil.moisture$Average.Soil.Moisture <- as.numeric(soil.moisture$Average.Soil.Moisture)
soil.moisture$Date <- as.factor(soil.moisture$Date)
is.na(soil.moisture)
moisture <- na.omit(soil.moisture)

#Create a new DF that groups by averages per plot
ASM.Plot <- moisture %>%
  group_by(Plot) %>%
  summarize(ASM=mean(Average.Soil.Moisture))

SM <- merge(moisture,ASM.Plot,by="Plot")
str(SM)
SM$Plot <- as.factor(SM$Plot)
#Data Visualizations
EL <- read.csv("Emerald Lake Plot Data.csv")

ggplot(EL, aes(x=Block, y=Soil.Moisture.2023)) +
  geom_point() + 
  geom_smooth(method=lm, se=TRUE)

#linear regression and ANOVA to to see if Block and Elevation effect SM
sm.lm <- lm(Soil.Moisture.2023~Block*Elevation, data = EL)
summary(sm.lm)
Anova(sm.lm)
