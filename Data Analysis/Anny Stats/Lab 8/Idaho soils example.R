#Idaho soils mixed model example
#AC 210226


#Set up---------------------------
#Set working directory, clean environment, etc.
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Anny Stats/Lab 8")

#load libraries
library(car)
library(lme4)#this is new
library(lmerTest)
library(performance)#this is new
library(see)#this is new
library(ggpubr)
library(emmeans)
library(patchwork)

#Import data and do some wrangling first------------------------
soil <- read.csv("Midwest labs soil analysis 190725.csv")

str(soil)
levels(soil$Sample_ID)#This variable contains the combinations of predictors


#Separate out our predictors from the "Sample.ID" variable (converts to multiple strings)
soil$Type<-substr(soil$Sample_ID,(nchar(as.character(soil$Sample_ID))-1),nchar(as.character(soil$Sample_ID)))
soil$Source<-substr(soil$Sample_ID,1,4)
soil$PlantID<-substr(soil$Sample_ID,1,nchar(as.character(soil$Sample_ID))-3)#Our random effect

#Now we're ready to actually do things
#Overall question, do soil sterilization and source affect soil chemistry response variables?
#We'll use a series of examples to see what happens when the random effect is included vs. not

#Example 1: Organic matter (OM)-----------------
#Prelim visualization
hist(soil$OM)
boxplot(soil$OM~soil$Type*soil$Source)

#Try to illustrate the influence of plantID
ggplot(soil, aes(x=Type, y=OM, group=PlantID,col=PlantID))+
  geom_point(alpha=0.5)+
  geom_line()+
  facet_wrap(vars(Source))

#Build a series of models to illustrate random effects
#1. No random effects
m.OM<-lm(OM~Type*Source,data=soil)
check_model(m.OM)
#2 Random intercept only
rI.OM<-lmer(OM~Type*Source+(1|PlantID),data=soil)
check_model(rI.OM)
#3 Random slope
rIS.OM<-lmer(OM~Type*Source+(0+Type|PlantID),data=soil)
#Does not work because we don't have a large enough sample size

#Compare outputs
summary(m.OM)
summary(rI.OM)
Anova(m.OM) 
Anova(rI.OM)

#Pairwise comparisions
emmeans(m.OM,pairwise~Source, adjust="Tukey")
emmeans(rI.OM,pairwise~Source, adjust="Tukey")

#Results figure
ggerrorplot(data=soil, x="Source", y="OM")


#Example 2: pH----------------------
#Prelim visualization
hist(soil$pH)
boxplot(soil$pH~soil$Type*soil$Source)

#Try to illustrate the influence of plantID
ggplot(soil, aes(x=Type, y=pH, group=PlantID,col=PlantID))+
  geom_point(alpha=0.5)+
  geom_line()+
  facet_wrap(vars(Source))

#Build a series of models to illustrate random effects
#1. No random effects
m.pH<-lm(pH~Type*Source,data=soil)
check_model(m.pH)
#2 Random intercept only
rI.pH<-lmer(pH~Type*Source+(1|PlantID),data=soil)
check_model(rI.pH)

#Compare outputs
summary(m.pH)
summary(rI.pH)
Anova(m.pH) 
Anova(rI.pH)

#Results figure
ggerrorplot(data=soil, x="Type", y="pH")

