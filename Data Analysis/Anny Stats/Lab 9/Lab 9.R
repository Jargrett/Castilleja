#Lab 9

#Set up---------------------------
#Set working directory, clean environment, etc.
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Data Analysis/Anny Stats/Lab 9")

#load libraries
library(car)
library(lme4)
library(lmerTest)
library(performance)
library(see)
library(ggpubr)
library(emmeans)
library(patchwork)
library(GGally)
library(factoextra)
library(tidyverse)

soil <- read.csv("Midwest labs soil analysis 190725.csv")

#Separate out our predictors from the "Sample.ID" variable (converts to multiple strings)
soil$Type<-substr(soil$Sample_ID,(nchar(as.character(soil$Sample_ID))-1),nchar(as.character(soil$Sample_ID)))
soil$Source<-substr(soil$Sample_ID,1,4)
soil$PlantID<-substr(soil$Sample_ID,1,nchar(as.character(soil$Sample_ID))-3)

#Question 3
ggpairs(soil, mapping = aes(color = Source), 
        columns = c( "OM","P1", "P2", 
                    "BICARB", "K", "Mg", "Ca", "pH", "Nitrate", "CEC"))+
  scale_colour_manual(values = c("darkorange","purple","cyan4", "blue")) +
  scale_fill_manual(values = c("darkorange","purple","cyan4", "blue"))

#Question 4
soil.pca <- prcomp (~ OM + P1 + P2 + BICARB + K + Mg + Ca + pH + Nitrate + CEC,
                    data=soil,
                    scale. = TRUE)
soil.pca
fviz_eig(soil.pca,addlabels = TRUE)

#Question 5
fviz_pca_biplot(soil.pca, label = "var",
                col.ind = soil$Source, palette = c("darkorange","purple","cyan4","blue"),
                col.var = "black", repel = TRUE,
                legend.title = "Source")
