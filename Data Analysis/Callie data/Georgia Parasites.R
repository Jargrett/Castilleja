#Georgia Parasites
#Initiated: 10/18/24

setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Callie data")

#load-in Packages
library(tidyverse)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plotting
library(car)#for regression analysis
library(emmeans)#post-hoc analysis
library(gridExtra)
library(webr)
library(formattable)
library(data.table)

#Load-in Data
georgia.parasites <- read.csv("Georgia Parasites - Master.csv")
georgia.parasites$taxa <- interaction(georgia.parasites$genus, georgia.parasites$species, sep = " ")
georgia.parasites <- georgia.parasites %>% relocate(taxa)
#Removing unnecessary coloumns
georgia.data <- georgia.parasites[ -c(32:36)]
#organize data strucutre, remove link columns, changing to factors
str(georgia.data)
georgia.data <- as.data.frame(unclass(georgia.data),stringsAsFactors=TRUE)
summary(georgia.data)
#Question 1:
# What is the current breakdown of status for Georgia's root hemiprasitic plants
# We only have data for 20 out of 40 Root HP:
# of which S1/S1? = 11, S2/S2? = 5, S3 = 1, SH = 1
filter.data <- georgia.data %>% 
  filter((parasitic_habit == "Root Hemiparasite")) %>% 
  filter((native == "Yes"))

summary(filter.data)

ggplot(filter.data, aes(x = ecosystem)) + 
  geom_bar() +
  geom_bar(aes(x = DNR_status_GA))








ecosystem.data <- plyr::rename(georgia.data, c("ecosystem_lump" = "Ecosystem",
                                               "habit" = "Habit",
                                               "DNR_status_lump" = "Status"))

PieDonut(ecosystem.data, aes(x = Habit, y = genus), r0=0.3, pieLabelSize = 4.7, donutLabelSize = 2.4,showRatioThreshold =.001, labelpositionThreshold=.01,
         ratioByGroup = TRUE)

PieDonut(ecosystem.data, aes(x = Ecosystem, y = habit), r0=0.3, pieLabelSize = 5, donutLabelSize = 2.4, explode=0,showRatioThreshold =.001, labelpositionThreshold=.01,
         ratioByGroup = TRUE)

PieDonut(ecosystem.data, aes(x = habit, y = NatureServ_status), r0=0.3, pieLabelSize = 5, donutLabelSize = 2.4, explode=0,showRatioThreshold =.001, labelpositionThreshold=.01,
         ratioByGroup = TRUE)

filter.data <- GP %>% 
  filter((parasitic_habit == "Root Hemiparasite")) %>% 
  filter((native == "Yes")) 

PieDonut(filter.data, aes(x = Status, y = DNR_monitor), r0=0.3, pieLabelSize = 5, donutLabelSize = 2.4,showRatioThreshold =.001, labelpositionThreshold=.01,
         ratioByGroup = TRUE)


write.csv(ecosystem.data, "/Users/jargrett/Desktop/Castilleja/Data Analysis/Callie data/GP.csv", row.names=FALSE)
GP <- read.csv("GP.csv")
GP <- as.data.frame(unclass(GP),stringsAsFactors=TRUE)

GP.filter <- GP %>% 
  filter((parasitic_habit == "Root Hemiparasite")) %>% 
  filter((native == "Yes")) %>% 
  filter((Status != "N/A"))

ggplot(GP.filter, aes(x = NatureServ_status, fill = DNR_monitor)) + 
  geom_bar()
  

