setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL Castilleja Observational Project")
#load in relevant packages
library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(performance)#this is new
library(see)#this is new
library(patchwork)
library(rstatix)
library(vegan)#for diversity analysis

#Diversity and Cover data file (combined in excel)
castilleja.cover <- read.csv("Processed Data/castilleja cover complete.csv")
castilleja.cover$castilleja[castilleja.cover$castilleja == "Control"] <- "Absent"
castilleja.cover$castilleja[castilleja.cover$castilleja == "Castilleja"] <- "Present"
castilleja.cover <- as.data.frame(unclass(castilleja.cover),stringsAsFactors=TRUE)

#Diversity calculations
#Import Datasets: We will be using the 23_24 Combined dataset for our analysis
case.cover <- read.csv("Raw Data/Case Cover 23_24 Combined - Cover.csv")
cali.cover <- read.csv("Raw Data/Cali Cover 23_24 Combined - Cover.csv")
cacr.cover <- read.csv("Raw Data/Cacr Cover 24 - Cover.csv")

#sperating the species matrix from the environmental data
#also will isolate bareground with the environmental data
cali.env <- subset(cali.cover, select=c(1:3,5:10))
case.env <- subset(case.cover, select=c(1:3,5:10))
cacr.env <- subset(cacr.cover, select=c(1:3,6:11))

#Isolating the species matrix, we also will remove Castilleja from the analysis to assess the background community
nocase.cover.matrix <- case.cover[ -c(1:10,26)]
nocali.cover.matrix <- cali.cover[ -c(1:10,25)]
nocacr.cover.matrix <- cacr.cover[ -c(1:11,30)]

#Isolating the species matric with Castilleja included we will use this to run the secondary analysis
case.cover.matrix <- case.cover[ -c(1:10)]
cali.cover.matrix <- cali.cover[ -c(1:10)]
cacr.cover.matrix <- cacr.cover[ -c(1:11)]
#--------------------Calculating Diversity Values--------------------#
#Using the vegan package we will calculate Shannon Diversity, Species Richness, and Pielou's evenness
#We are using our cover data for this analysis, we will use the conservative no Castilleja matrix 

#cali
cali.cover.div <- diversity(nocali.cover.matrix, index = "shannon")
cali.cover.rich <- specnumber(nocali.cover.matrix)
cali.cover.even <- diversity(nocali.cover.matrix, index = "shannon") / log(specnumber(nocali.cover.matrix)) 
#case
case.cover.div <- diversity(nocase.cover.matrix, index = "shannon")
case.cover.rich <- specnumber(nocase.cover.matrix)
case.cover.even <- diversity(nocase.cover.matrix, index = "shannon") / log(specnumber(nocase.cover.matrix))
#cacr
cacr.cover.div <- diversity(nocacr.cover.matrix, index = "shannon")
cacr.cover.rich <- specnumber(nocacr.cover.matrix)
cacr.cover.even <- diversity(nocacr.cover.matrix, index = "shannon") / log(specnumber(nocacr.cover.matrix))

#Combine calculated values with our environmental data
cali.cover.diversity <- cbind(cali.env,cali.cover.div,cali.cover.rich,cali.cover.even)
case.cover.diversity <- cbind(case.env,case.cover.div,case.cover.rich,case.cover.even)
cacr.cover.diversity <- cbind(cacr.env,cacr.cover.div,cacr.cover.rich,cacr.cover.even)
#renaming columns
case.cover.diversity <- plyr::rename(case.cover.diversity, c("case.cover.div" = "div",
                                                             "case.cover.even" = "even",
                                                             "case.cover.rich" = "rich"))

cali.cover.diversity <- plyr::rename(cali.cover.diversity, c("cali.cover.div" = "div",
                                                             "cali.cover.even" = "even",
                                                             "cali.cover.rich" = "rich"))

cacr.cover.diversity <- plyr::rename(cacr.cover.diversity, c("cacr.cover.div" = "div",
                                                             "cacr.cover.even" = "even",
                                                             "cacr.cover.rich" = "rich"))
#creating a new csv 
write.csv(cali.cover.diversity, "Processed Data/Combined linariifolia Cover Diversity.csv", row.names=FALSE)
write.csv(case.cover.diversity, "Processed Data/Combined septentrionalis Cover Diversity.csv", row.names=FALSE)
write.csv(cacr.cover.diversity, "Processed Data/Chromosa Cover Diversity.csv", row.names=FALSE)

