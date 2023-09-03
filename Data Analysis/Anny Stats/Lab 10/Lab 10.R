
#Set up---------------------------
#Set working directory, clean environment, etc.
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Data Analysis/Anny Stats/Lab 10")
#Load packages
library(vegan)#this is new
library(ggplot2)
library(ggrepel)

#Load data--------------------------------
fire <-read.csv("sagefire_abund_matrix.csv", row.names=1)
#Matrix has 70 samples (rows) with abundances of 35 mite species
#Note that the "site names" have to be input as row names. Cannot be separate column

fire.info<-read.csv("sagefire_info.csv")

fire.info$Firefreq<-as.factor(fire.info$FireFreq)

set.seed(20)

fire.dist<-vegdist(fire, method="bray")

fire.nmds<-metaMDS(fire.dist,distance="bray", k=2,try=100)
fire.nmds
stressplot(fire.nmds)

#PERMANOVA
adonis(fire.dist~Firefreq, data = fire.info, permutations=999)

#plot
nmds.scores <- as.data.frame(vegan::scores(fire.nmds))
fire.info<-cbind.data.frame(fire.info,nmds.scores) 
fire.info$site <- row.names(fire)

ggplot(fire.info, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=Firefreq, shape=Firefreq)) +
  geom_text_repel(aes(NMDS1, NMDS2, label = site), size = 2.5)

