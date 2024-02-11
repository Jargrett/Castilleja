# Deer Creek Site 1 Cover Analysis
# Initiated: 7/16/23
# Completed: TBD

####################
#   Fixed Effects  # - Treatment (Presence/Absence)
#Response Variables# - Cover
#  Random Effects  # - Plot
#     key Terms    # - 
####################

#####################
#                   # - 
#      Question     # - 
#                   # - 
##################### - 

# Set working directory (Workspace)
# This can be done manually or through session -> set working directory
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Data Analysis/RMBL/Makena Data")

# Load-in packages
library(car)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(lme4)
# I will explain these if needed

#Load data
dc1.cover <- read.csv("dc site 1 cover.csv")

# It is always a good idea to view the dataset to see if its correct
View(dc1.cover)

#---------------------NOW WE BEGIN!!!!!---------------------

# Checking data structure
# this allows you to see what is coded as 
str(dc1.cover)

# change these chars/integers to factors so that we can use them in analysis
dc1.cover$Treatment <- as.factor(dc1.cover$Treatment)
dc1.cover$Plot <- as.factor(dc1.cover$Plot)

#More data viewing (Look at how zero inflated the data is, this is ok!!)
table(dc1.cover$Castilleja.linariifolia)
table(dc1.cover$Bare.ground)

# We can use tidyverse (package) to shift our dataframe
# We will combine all species into one column and cover into another while maintaining the structure of or original dataset for the rest
dc1.cover.long <- pivot_longer(dc1.cover, col = 5:42, names_to = "species", values_to = "cover")

# Just to show you, that our data is very 0 inflated so we cant really see anything from this alone
boxplot(dc1.cover.long$cover~dc1.cover.long$Treatment)

# Removing 0 values by first changing them to NA
dc1.cover.long[dc1.cover.long==0] <- NA
dc1.cover.zero<-dc1.cover.long[complete.cases(dc1.cover.long),]

# lets look again
boxplot(dc1.cover.zero$cover~dc1.cover.zero$Treatment)#looks a bit better
str(dc1.cover.zero)

# Changing species to a factor
dc1.cover.zero$species <- as.factor(dc1.cover.zero$species)
str(dc1.cover.zero)

#now lets create a model that looks at the relationship between total cover and treatment
cm<-lm(log(cover)~Treatment, data=dc1.cover.zero)
hist(cm$residuals)
qqPlot(cm)#much better
plot(cm)#seems okay

# We will now run and analysis of variance(ANOVA) to assess the effect or treatment
Anova(cm)#Not significant 0.8213
#Bummer for now but 

#--------------------Diversity Time!!--------------------#
# For this we will go back to dc1.cover as it is in the correct format
# this is why its always good to save your alterations into a new data file rather than altering an old one
#lets look at it again

View(dc1.cover) # looks almost ready to go minus a few things

dc1.diversity <- dc1.cover[ -c(1,2,4) ]
dc1.cover.diversity$Plot <- as.numeric(dc1.cover.diversity$Plot)
str(dc1.cover.diversity)
diversity(dc1.cover.diversity,index = "shannon")
