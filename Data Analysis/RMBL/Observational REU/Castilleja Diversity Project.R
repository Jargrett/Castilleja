#Castilleja Diversity Project
#Initiated: 8/27/24

#Set the working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU")

#load in relevant packages
library(tidyverse)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring

#Import Datasets: We will be using the 23_24 Combined dataset for our analysis
case.cover <- read.csv("Case Cover 23_24 Combined - Cover.csv")
cali.cover <- read.csv("Cali Cover 23_24 Combined - Cover.csv")

#sperating the species matrix from the environmental data
#also will isolate bareground with the environmental data
cali.env <- subset(cali.cover, select=c(1:3,5:10))
case.env <- subset(case.cover, select=c(1:3,5:10))

#Isolating the species matrix, we also will remove Castilleja from the analysis to assess the background community
nocase.cover.matrix <- case.cover[ -c(1:10,26)]
nocali.cover.matrix <- cali.cover[ -c(1:10,25)]
#Isolating the species matric with Castilleja included we will use this to run the secondary analysis
case.cover.matrix <- case.cover[ -c(1:10)]
cali.cover.matrix <- cali.cover[ -c(1:10)]

#--------------------Calculating Diversity Values--------------------#
library(vegan)#for diversity analysis
#Using the vegan package we will calcuate Shannon Diversity, Species Richness, and Pielou's evenness
#We are using our cover data for this analysis, we will use the conservative no Castilleja matrix 

#cali
cali.cover.div <- diversity(nocali.cover.matrix, index = "shannon")
cali.cover.rich <- specnumber(nocali.cover.matrix)
cali.cover.even <- diversity(nocali.cover.matrix, index = "shannon") / log(specnumber(cali.cover.matrix)) 
#case
case.cover.div <- diversity(nocase.cover.matrix, index = "shannon")
case.cover.rich <- specnumber(nocase.cover.matrix)
case.cover.even <- diversity(nocase.cover.matrix, index = "shannon") / log(specnumber(case.cover.matrix))

#Combine calculated values with our environmental data
cali.cover.diversity <- cbind(cali.env,cali.cover.div,cali.cover.rich,cali.cover.even)
case.cover.diversity <- cbind(case.env,case.cover.div,case.cover.rich,case.cover.even)

#renaming columns
case.cover.diversity <- plyr::rename(case.cover.diversity, c("case.cover.div" = "div",
                                                             "case.cover.even" = "even",
                                                             "case.cover.rich" = "rich"))

cali.cover.diversity <- plyr::rename(cali.cover.diversity, c("cali.cover.div" = "div",
                                                             "cali.cover.even" = "even",
                                                             "cali.cover.rich" = "rich"))

#creating a new csv for later use and sharing
write.csv(cali.cover.diversity, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined linariifolia Cover Diversity.csv", row.names=FALSE)
cali.cover.diversity <- read.csv("Combined linariifolia Cover Diversity.csv")
write.csv(case.cover.diversity, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/Observational REU/Combined septentrionalis Cover Diversity.csv", row.names=FALSE)
case.cover.diversity <- read.csv("Combined septentrionalis Cover Diversity.csv")

#--------------------Diversity Analysis--------------------#
#We will do analysis by species starting with septentrionalis first

#Lets first check the data structure
str(case.cover.diversity)
#We will want to change some of our integers and characters to factors
#this will allow us to group and compare by these factor groupings

case.cover.diversity$site <- as.factor(case.cover.diversity$site)
case.cover.diversity$year <- as.factor(case.cover.diversity$year)
case.cover.diversity$castilleja <- as.factor(case.cover.diversity$castilleja)
case.cover.diversity$pair <- as.factor(case.cover.diversity$pair)
case.cover.diversity$plot <- as.factor(case.cover.diversity$plot)

#Now we can begin by running a linear mixed effect model
#Fixed effects = castilleja and site
#Random effects = year and pair

#Shannon Diverisity
case.div <- lmer(div ~ castilleja*site + (1|year) + (1|pair), data = case.cover.diversity)
summary(case.div)
Anova(case.div) #Castilleja p = 0.02415, Site p =  < 2e-16
emmip(case.div, castilleja ~ site)
emmeans(case.div, pairwise ~ castilleja|site) #Seems that difference is driven by EL

#Species Richness
case.rich <- lmer(rich ~ castilleja*site + (1|year) + (1|pair), data = case.cover.diversity)
summary(case.rich)
Anova(case.rich) #Castilleja p = 0.002825, Site p = < 2.2e-16
emmip(case.rich, castilleja ~ site)
emmeans(case.rich, pairwise ~ castilleja|site) #Seems that difference is again driven by EL

#Species Richness
case.even <- lmer(even ~ castilleja*site + (1|year) + (1|pair), data = case.cover.diversity)
summary(case.even)
Anova(case.even) #Castilleja p = 0.0001496, Site p = < 2.2e-16, insignificant interaction p = 0.0878680
emmip(case.even, castilleja ~ site)
emmeans(case.even, pairwise ~ castilleja|site) #Seems that difference is again driven by Avery control is more even


#Now we will do the same for linariifolia, code should be more or less the same
cali.cover.diversity$site <- as.factor(cali.cover.diversity$site)
cali.cover.diversity$year <- as.factor(cali.cover.diversity$year)
cali.cover.diversity$castilleja <- as.factor(cali.cover.diversity$castilleja)
cali.cover.diversity$pair <- as.factor(cali.cover.diversity$pair)
cali.cover.diversity$plot <- as.factor(cali.cover.diversity$plot)

#Now we can begin by running a linear mixed effect model
#Fixed effects = castilleja and site
#Random effects = year and pair

#Shannon Diverisity
cali.div <- lmer(div ~ castilleja*site + (1|year) + (1|pair), data = cali.cover.diversity)
summary(cali.div)
Anova(cali.div) #Castilleja p = 0.006585, Castilleja by site interaction p = 0.018733, insignificant site p = 0.052747
emmip(cali.div, castilleja ~ site) # looks like magnitudes are different between sites
emmeans(cali.div, pairwise ~ castilleja|site) #Seems that difference is driven by DC2

#Species Richness
cali.rich <- lmer(rich ~ castilleja*site + (1|year) + (1|pair), data = cali.cover.diversity)
summary(cali.rich)
Anova(cali.rich)#Castilleja p = 00.006421, Castilleja by site interaction p = 0.007023
emmip(cali.rich, castilleja ~ site)# magnitude differences as well as Johnson inversing the relationship
emmeans(cali.rich, pairwise ~ castilleja|site) #Seems that difference is again driven by DC2

#Species Richness
cali.even <- lmer(even ~ castilleja*site + (1|year) + (1|pair), data = cali.cover.diversity)
summary(cali.even)
Anova(cali.even) #Site p = 0.003408
emmip(cali.even, castilleja ~ site)# lot going on here, looks like no consistant trend between sites
emmeans(cali.even, pairwise ~ castilleja|site)#DC1 difference is sig p = 0.0223 but the retalionship is higher in control similar to Case data


