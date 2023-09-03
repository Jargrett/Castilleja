# Emerald Lake Pre-Treatment plant cover
# Initiated: 7/25/23
# Completed: TBD

#####################
#                   # - Cover analysis of paired plots
#      Question     # - Total Castilleja cover
#                   # - NN analysis and richness
##################### - 

# Set working directory
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Data Analysis/RMBL/2023 Argrett")

# Load-in packages
library(tidyverse) # for data working
library(vegan) # for diversity
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car)
library(ggpubr)
library(emmeans) # for comparison of means
library(rstatix) # for comparison of means
library(labdsv)

# Loading in data
plants <- read.csv("Pre Treatment Plant Data.csv")
#removing NA columns
##plants <- pretreatment %>% select_if(~ !any(is.na(.)))


#-----------------Checking data structure-----------------#
str(plants)
plants$Code <- as.factor(plants$Code)
plants$Functional.Group <- as.factor(plants$Functional.Group)
plants$Block <- as.factor(plants$Block)
plants$Removal <- as.factor(plants$Removal)
plants$Life.History <- as.factor(plants$Life.History)

#------------Diversity Analysis------------#
# here we will be working with the count column
# we will need to conververt data to a matrix format
plant.c <- subset(plants, select = c('Plot','Code','Count'))
plant.count<- matrify(plant.c)

# Calculating the species richness for plots
pre.rich <- specnumber(plant.count)

# Calculating Shannon diversity for plots
pre.div <- diversity(plant.count, index = "shannon")

#combined data set with Plot Data File and calculated values
plot.data <- read.csv("Plot Data.csv") #importing metadata
p.data <- plot.data %>% select_if(~ !any(is.na(.))) #removing NA Columns
plant.div <- cbind(p.data, pre.rich, pre.div) #final datset

#Run some models
pre.rich.glm <- glm(pre.rich ~ Block*Elevation, family="poisson", data = plant.div)
hist(pre.rich.glm$residuals) #pretty normal
plot(pre.rich.glm)#seems okay

# We can get a summary of the model:
summary(pre.rich.glm) 

pre.rich.plot <- ggplot(plant.div, aes(x = Elevation, y = pre.rich)) +
  geom_point() +
  labs(x = "Elevation", y = "Species Richness") +
  theme_bw()
pre.rich.plot


#---------------------Nearest Neighbor Analysis---------------------#
# Count data for nearest neighbor needs to be cleaned and processed
# End result is a table with |Species|Neighbor Count|Cover|Haustoria Found|

#subset the dataset
pre.neighbor <- subset(plants,select=c('Plot','Block','Species','Code','Functional.Group','Nearest.Neighbor'))
plant.neighbor <- filter(pre.neighbor, Nearest.Neighbor > 0) # filter 0 rows
table(plant.neighbor$Functional.Group)


# visualize
plot.neighbor<-ggplot(data=plant.neighbor, aes(x= Functional.Group, y= Nearest.Neighbor)) +
  geom_bar(stat="identity")

plot.neighbor


#---------------------Cover Analysis---------------------#
# What are the top 5 most dominant plants by cover?
# CASE(123.6) CHAN(268.8) Elymus(260.3) HEMA(173.6) HEQU(173.6)
# LIPO(275.5) MESP(104.3) POGR(112.4) RIMO(199) HYHO(172.5)

table(plants$Disturbance.Type)
table(plants$Block)
# how is cover affected by elevation and disturbance?
# to subset =subset(myData, Code =="CASE")
m <- lm(percent.cover~Block*Elevation, data = plants)
summary(m)

#Elevation = positive effect on cover of RIMO(.033), Elymus(.009)
#Elevation = negative effect on cover of CHAN(.004)

plot.species.cover<-ggplot(plants[plants$Code %in% "HEQU" ], aes(x=Elevation, y=percent.cover)) +
  geom_point() +
  geom_smooth(method=lm)
plot.species.cover

plot.species.cover<-ggplot(plants, aes(x=Block, y=percent.cover)) +
  geom_boxplot() +
  geom_smooth(method=lm)
plot.species.cover

