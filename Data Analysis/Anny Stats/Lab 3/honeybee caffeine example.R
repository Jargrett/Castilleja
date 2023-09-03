#One-way ANOVA example: honeybees
#AC 210112

#Do your housekeeping----------------------
#Set working directory
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Anny Stats/Lab 3")
#Make sure your environment is clean

#Load libraries for this exercise---------
library(car)
library(ggplot2)
library(ggpubr)

#Load data and check---------------------------------
bees<-read.csv("honeybee.csv")

#Check data structure
str(bees)
#Change ppmCaffeine to factor
bees$ppmCaffeine <- as.factor(bees$ppmCaffeine)

#Preliminar data visualization-----------
#Make histogram of response variable (turn in)
hist(bees$consumptionDifferenceFromControl)
#Make boxplots of response~treatment levels (turn in)
boxplot(bees$consumptionDifferenceFromControl~bees$ppmCaffeine)
#Define linear model and check assumptions----------
bm1<-lm(consumptionDifferenceFromControl~ppmCaffeine, data=bees)

#Check assumptions:
#Make qqplot of residuals (turn in)
qqPlot(bm1)
#Check homogeneity of variance (turn in)
plot(bm1)

#Examine results of linear model-------------------
#1. summary()
summary(bm1)
#How to figure out what "Intercept" represents
bees$ppmCaffeine#R automatically chooses the first 1 (usually lowest number or alphabetically) as baseline
#the model matrix shows how the effects are defined/calculated
model.matrix(bm1)

#2. anova()
anova(bm1)

#3. summary(aov())
summary(aov(consumptionDifferenceFromControl~ppmCaffeine, data=bees))

#4. Anova() in car package
Anova(bm1)

#Post-hoc comparisons of means (Tukey-HSD)----------
#Note TukeyHSD() only works on aov objects
baov <- (aov(consumptionDifferenceFromControl~ppmCaffeine, data=bees))
TukeyHSD(baov)



#Make a mean and standard error plot ----------------
Beeplot <- ggerrorplot(bees, x = "ppmCaffeine" , y = "consumptionDifferenceFromControl", 
            desc_stat = "mean_se", add.params = list(color = "darkgray"),
            ylab  = "Consumption different from Control", xlab = "Caffine concentration (ppm)" )

Beeplot

#Suggest using ggpubr, but use whatever method works for you