#Set the working directory
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Callie Data/CURO")

#load in relevant packages
library(plyr)#for data wrangling and restructuring
library(tidyverse)#for data wrangling and restructuring
library(statmod)

biomass <- read.csv("Callie Data - Biomass.csv")
biomass <- as.data.frame(unclass(biomass),stringsAsFactors=TRUE)


biomass.lm <- lm(SOAL_biomass ~ treatment*parasite, data = biomass)
summary(biomass.lm)
Anova(biomass.lm)
emmeans(biomass.lm, pairwise ~ treatment|parasite)

biomass.lm <- lm(AGPU_biomass ~ treatment, data = biomass)
summary(biomass.lm)
Anova(biomass.lm)
emmeans(biomass.lm, pairwise ~ treatment|parasite)


#this is scratch
totCmod<-lmer(totC.bulk~Tree*Burn+propClay+Depth+pH+(1|Site), data=op)
Cemmeans<-emmeans(totCmod, specs=pairwise~Tree|Burn)$emmeans
C_means<-summary(Cemmeans)$emmean
C_SE<-summary(Cemmeans)$SE
dfC<-data.frame(
  Tree=factor(c("Oak","Pine","Oak","Pine")),
  Burn=factor(c("Unburned","Unburned","Burned","Burned")),
  mean=C_means,
  lower=C_means-C_SE,
  upper=C_means+C_SE
)
c_lett<-cld(Cemmeans, decreasing=T, Letters=letters, sort=T)
ggplot(data=dfC)+geom_pointrange(aes(x=Tree,y=mean, ymin=lower,ymax=upper))+
  facet_wrap(vars(Burn))+
  geom_text(aes(label = c_lett$.group, x=c(1.1,2.1,1.1,2.1), y=c(2,2,2,2)))