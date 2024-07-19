#Castilleja Diveristy Analysis

#Set the working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU")

#load in packages
library(tidyverse)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(lme4)#for modeling linear mixed effect models
library(nlme)#alternative for modeling linear mixed effect models
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plotting
library(car)#for regression analysis
library(vegan)#for diversity analysis
library(emmeans)#post-hoc analysis

#Load in 2023 + 2024 datasets (For diversity analysis we will use species counts)
CALI.24 <- read.csv("Cali 2024 Cover.csv")
CALI.24.Johnson <- read.csv("Cali 2024 Johnson Hill Cover.csv")
CALI.23 <- read.csv("Cali 2023 Cover.csv")
CASE.23 <- read.csv("Case 2023 Cover.csv")
CASE.24.Avery <- read.csv("Case 2024 Avery Cover.csv")

#merging dataframes by pivot longer then rbind
cali <- rbind.fill(CALI.23,CALI.24,CALI.24.Johnson)
case <- rbind.fill(CASE.23,CASE.24.Avery)
combined <- rbind.fill(cali,case)

#changing NA to 0
cali[is.na(cali)] <- 0
case[is.na(case)] <- 0
combined[is.na(combined)] <- 0


# str(CASE.24.Avery)
# cali.24.long<- pivot_longer(CALI.24, cols = Achratherum.lettermanii:Wyethia.amplexicaulis,
#                          names_to = "species",
#                          values_to = "count")
# cali.24.johnson.long<- pivot_longer(CALI.24.Johnson, cols = Achnatherum.sp.:Viola.praemorsa,
#                             names_to = "species",
#                             values_to = "count")
# cali.23.long<- pivot_longer(CALI.23, cols = Agastache.urticifolia:Wyethia.amplexicaulis,
#                                     names_to = "species",
#                                     values_to = "count")
# case.23.long<- pivot_longer(CASE.23, cols = Agoseris.glauca.leaves:Viola.adunca,
#                             names_to = "species",
#                             values_to = "count")
# case.24.Avery.long<- pivot_longer(CASE.24.Avery, cols = Achillea.millefolium:Viola.adunca,
#                             names_to = "species",
#                             values_to = "count")
# 
# combined.counts <- rbind(cali.23.long,cali.24.long, DF3)

