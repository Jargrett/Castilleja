#Final Project ROUGH CODE

########################
#!Use Lots of Comments!#
########################

#setwd
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Data Analysis/Anny Stats/Final Project")

#Load Packages
library(car)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(emmeans)
library(lme4)
library(nlme)
library(MuMIn)
library(multcomp)


#Data Link
## https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-sev.186.208430 

#Load data
SEV <- read.csv("sev186_NPP_fertilizer_biomass.csv")

####################
#   Fixed Effects  # - treatment, photopath, year, lifeHistory
#Response Variables# - cover, Biomass.BM
#  Random Effects  # - season, season.precip, plot, quad
#     key Terms    # - karetz (species code), Biomass.BM (best fit)
####################

#####################
#                   # - impact of N treatment on biomass of c4 vs c3 species over time?
#      Question     # - impact of N treatment on biomass of annuals vs perrenials
#                   # - both within functional group photopathways as well as overall
##################### - how does photopath and functional group cause biomass/cover shift

#Checking structure and data cleaning
str(SEV)
SEV$year <- as.factor(SEV$year)
# We have c3,c4,CAM, and one mixed c3/c4
table(SEV$PhotoPath)
# We have annuals, perennials, and a mix
table(SEV$LifeHistory)
# We have both a control and a fertilizer treatment
table(SEV$treatment)
#experiment runs from 2004 to 2021 (18 yrs)
table(SEV$year) # why do you final years have less counts then before??
#20 total plots with 4 nested quadrats in a plot
table(SEV$plot)
table (SEV$quad)
#data is 0 skewed due to it being biomass data
boxplot(SEV$biomass.BM~SEV$treatment)

#Subsetting data and further data modification
#SEVD <- select(SEV, c(1,2,8,10:12,15:19,21:22))

SEVD <- SEV %>%
  dplyr::select(year, season, season.precip, plot, quad, treatment, kartez, family, LifeHistory, PhotoPath, FunctionalGroup, cover, biomass.BM) %>%
  filter(PhotoPath != "CAM" & PhotoPath != "C3-C4") 


unique(SEVD$PhotoPath)


SEVD$Stdyyear <- as.numeric(SEVD$year)
str(SEVD$Stdyyear)
SEVD0 <- filter(SEVD, PhotoPath != "CAM" & PhotoPath != "C3-C4") #removes CAM and C3/C4 Plants
SEVD1 <- filter(SEVD0, LifeHistory != "ann/perenn") #removes annual/perennial

# Let's take out last two years
SEVD2 <- filter(SEVD1, Stdyyear != 18)
SEVD2 <- filter(SEVD2, Stdyyear != 17)
table(SEVD2$Stdyyear)

# convert treatment and photo path to a factor
SEVD2$treatment <- as.factor(SEVD2$treatment)
SEVD2$PhotoPath <- as.factor(SEVD2$PhotoPath)


#log transform for cover and biomass

which(SEVD2$biomass.BM==0)
Bmin.nz <- min(SEVD2$biomass.BM[which(SEVD2$biomass.BM>0)]) #min non-zero

hist(SEVD2$log.biomass.BM) #looks much better

#I want to find the average biomass by year
table(SEVD2$PhotoPath, SEVD2$Stdyyear, SEVD2$treatment) #breakdown of C3/C4 counts per year
table(SEVD2$LifeHistory, SEVD2$Stdyyear, SEVD2$treatment) #breakdown of annual/perennial

#creating a new dataset with sum biomass by year
SEV.sum <- aggregate(biomass.BM~Stdyyear*plot*quad*PhotoPath*LifeHistory*season.precip*treatment, data = SEVD2, FUN = sum)
SEV.sum$treatment <- gsub("C", "Control", SEV.sum$treatment)
SEV.sum$treatment <- gsub("F", "Treatment", SEV.sum$treatment)
SEV.sum$log.biomass.BM <- log(SEV.sum$biomass.BM + 0.00066)
SEVD2$log.biomass.BM <- log(SEVD2$biomass.BM + 0.00066)
#More Data checks
plot(SEV.sum$Stdyyear, SEV.sum$log.biomass.BM) # not many trends that I can see throughout
plot(SEV.sum$Stdyyear, SEV.sum$season.precip)


#_______________________________Modeling_______________________________#
#---Time Series analysis of Biomass by treatment
bio_mod<-lme(log.biomass.BM~treatment*Stdyyear, random=~1|plot, correlation=corAR1(form=~1|plot), data=SEV.sum)
summary(bio_mod)
Anova(bio_mod)
##emmeans(bio_mod, pairwise ~treatment)
#---Time Series analysis of Biomass by photopah Q1
Pbio_mod<-lme(log.biomass.BM~treatment*PhotoPath*Stdyyear, random=~1|plot, correlation=corAR1(form=~1|plot), data=SEV.sum)
summary(Pbio_mod)
Anova(Pbio_mod)
##emmeans(Pbio_mod, pairwise ~treatment*PhotoPath)

#---Time Series analysis of Biomass by lifehistory Q2
Lbio_mod<-lme(log.biomass.BM~treatment*LifeHistory*Stdyyear, random=~1|plot, correlation=corAR1(form=~1|plot), data=SEV.sum)
summary(Lbio_mod)
Anova(Lbio_mod)
##emmeans(Lbio_mod, pairwise ~treatment*LifeHistory)
#---Time Series analysis of Biomass by photopth and lifehistory
Tbio_mod<-lme(log.biomass.BM~treatment*PhotoPath*LifeHistory*Stdyyear, random=~1|plot, correlation=corAR1(form=~1|plot), data=SEV.sum)
summary(Tbio_mod)
Anova(Tbio_mod)
##emmeans(Tbio_mod, pairwise ~treatment*PhotoPath*LifeHistory*Stdyyear)

#___________________________________Plots______________________________________#
#timeseries plot of treatments
timeseries <- ggplot()+
  geom_line(data=SEV.sum, 
          aes(x=Stdyyear,y=log.biomass.BM, group=factor(PhotoPath), linetype=PhotoPath),
          stat = "summary", fun = sum, size = 1) +
  facet_grid(LifeHistory ~ treatment) +
  labs(y = expression(paste("Log Aboveground biomass (g)")), x = expression(paste("Study Year"))) +
  theme_pubclean()
 
timeseries

#changes in biomass by treatment throughout entire study(not by year)

#mean seasonal precipitation by year 
#we can see that there is high variablilty in precip during the fall

plot2 <- ggplot()+
  geom_line(data = SEV.sum,
                 aes(x = Stdyyear, y = season.precip, group=factor(season), color=season),
                 stat = "summary", fun = mean, size = 0.8) +
  theme_pubclean()
plot2

#above ground biomass and life history over study year
#annuals have more fluctuation and variance
plot3 <- ggplot()+
  geom_line(data=SEVD2, 
            aes(x=Stdyyear,y=biomass.BM, group=factor(LifeHistory), color=LifeHistory),
            stat = "summary", fun = mean, size = .8) +
  facet_grid(treatment ~.) +
  labs(y = expression(paste("Average Aboveground biomass (g)")), x = expression(paste("Study Year"))) +
  theme_pubclean()

plot3

        