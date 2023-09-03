#Linear regression example: dragons
#AC 210110

#Do your housekeeping----------------------
#Set working directory
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Anny Stats/Lab 3")
#Make sure your environment is clean

#Load libraries for this exercise---------
library(car)
library(ggplot2)

#Load data---------------------------------
dragon <- read.csv("C:/Users/jargr/Dropbox/PC/Desktop/Anny Stats/Lab 3/dragons.csv")
View(dragons)

#Examine data variables--------------------

#In this example, we are trying to understand what predicts dragon weight
#We have two predictor variables, one continuous (height), one categorical (pattern)
#(assume that height is measured without error)

#Set "striped" as baseline manually
dragon$pattern = relevel(dragon$pattern, ref="striped")

#Preliminary visualization-----------------
#Histograms: what do the variable distributions look like
hist(dragon$weight)
#Box plot: understand general pattern, look out for ourliers and heteroscedasticity
boxplot(dragon$weight~dragon$pattern)
#X-Y plot: visualize relationship, make sure linear is reasonable
plot(x = dragon$height, y = dragon$weight)
#Conduct multiple regression---------------
#Some of your past stats classes probably call this an ANCOVA
#Define model (first part defines relationship between variables)
dm1 <- lm(weight~height+pattern, data = dragon )
#Check for violations of assumptions (normality and homoscedasticity)
plot(dm1)
qqplot(dm1)
hist(dm1)
summary(dm1)
#Take a look at the results (back to ppt)

#Visualize results----------------------
#Basic scatter plot
ggplot(data = dragon, aes(x=height, y=weight, color=pattern))+
  geom_point()+
  theme_classic()

#The safest way to add your fit lines is to use predicted values
#Do not rely on "reg.line" or "geom_smooth with lm"
#Add predicted values to data (fitted fits our model to data)
dragon$pred.weight<-fitted(dm1)

#Try plotting again, add lines using predicted values
ggplot(data = dragon, aes(x=height, y=weight, color=pattern))+
  geom_point()+
  geom_line(aes(x=height, y=pred.weight, color=pattern))+ #make sure to specify that lines are plotted with predicted weight
  theme_classic()

#How to make this plot better?
#Hint: units? colors? point sizes and shapes?

