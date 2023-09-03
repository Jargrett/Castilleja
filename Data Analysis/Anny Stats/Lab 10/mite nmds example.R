#NMDS in class example code 
#AC 210327

#Set up---------------------------
#Set working directory, clean environment, etc.
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Data Analysis/Anny Stats/Lab 10")
#Load packages
library(vegan)#this is new
library(ggplot2)
library(ggrepel)

#Load data--------------------------------
mite<-read.csv("mite_abund_matrix.csv", row.names=1)
#Matrix has 70 samples (rows) with abundances of 35 mite species
#Note that the "site names" have to be input as row names. Cannot be separate column

mite.info<-read.csv("mite_explain_var.csv")

#Conduct NMDS------------------------------
set.seed(20)#This sets the random start seed so that we are guaranteed to all get the same outputs

#First calculate distance matrix
mite.dist<-vegdist(mite, method="bray")


#Run NMDS on distance matrix
mite.nmds<-metaMDS(mite.dist,distance="bray", #use bray-curtis distance
                  k=2, #2 dimensions
                  try=100 #force it to try 100 different times (default is 20, for publication I recommend 500)
                  )
mite.nmds #Gives stress value of best solution, aim for <0.2

#Check the fit
stressplot(mite.nmds) #Also called a "Shepard diagram"
#We want this plot to show a monotonic relationship

#Make ordination plot--------------------------
#Using basic vegan graphics
#Site labels only
ordiplot(mite.nmds, type="text", display="sites")
#Points only
ordiplot(mite.nmds, type="points", display="sites")

#Can output ordination coordinates and use ggplot for nicer graphics
nmds.scores <- as.data.frame(vegan::scores(mite.nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
#Add to mite.info dataframe
mite.info<-cbind.data.frame(mite.info,nmds.scores) 
mite.info$site <- row.names(mite)

#Build NMDS figure in ggplot
ggplot(mite.info, aes(NMDS1, NMDS2)) +
  geom_point() +
  geom_text_repel(aes(NMDS1, NMDS2, label = site), size = 3)#avoids text overlapping

#Diagnostic detour: what happens if stress is too high?-----
#Usually the solution is increase k
#Compare previous example using k=3
mite.nmds3<-metaMDS(mite.dist,distance="bray", k=3,try=100)
mite.nmds3 #stress is lower than k=2, but that means you'd have to make a 3D plot to visualize

#Here's a nifty function that fits multiple k's to help visualize change in stress
#Also called a Scree plot
#define function to perform a NMDS for 1-6 dimensions and plots the nr of dimensions vs the stress
NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 6), replicate(6, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 6),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS Scree plot")
  for (i in 1:6) {
    points(rep(i + 1,6),replicate(6, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

#Try it on our data
# Use the function that we just defined to choose the optimal nr of dimensions
NMDS.scree(mite.dist)

#Example: Color NMDS plot by different shrub cover----------------------
#Revise plotting code from above
ggplot(mite.info, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=Shrub, shape=Shrub)) +
  geom_text_repel(aes(NMDS1, NMDS2, label = site), size = 3)#avoids text overlapping

#PERMANOVA to test whether composition differs by shrub cover or topography---------
adonis(mite.dist~Shrub*Topo, data = mite.info, permutations=999)
#The accuracy/significant figures of p value depends on permutations

#Make nice NMDS figure to reflect this PERMANOVA analysis------------------
ggplot(mite.info, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=Shrub, shape=Topo)) 

