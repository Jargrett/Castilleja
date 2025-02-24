#Penguins PCA example
#AC 210319

#Set up---------------------------
#Set working directory, clean environment, etc.
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Anny Stats/Lab 9")
#Load packages
library(GGally)#this is new
library(factoextra)#this is new
library(tidyverse)


#Import data--------------------
penguins <- read.csv("penguins.csv")

#remove rows if NAs occur in any of the variables we care about
penguins1 <- penguins %>% filter_at(vars(flipper_length_mm, body_mass_g,bill_length_mm, bill_depth_mm),
                                    all_vars(!is.na(.)))

#Examine correlations among the body shape response variables------
#ggpairs plots is a good way to get lots of info (GGally package)
ggpairs(penguins1, mapping = aes(color = species), 
        columns = c("flipper_length_mm", "body_mass_g", 
                    "bill_length_mm", "bill_depth_mm"))+
  scale_colour_manual(values = c("darkorange","purple","cyan4")) +
  scale_fill_manual(values = c("darkorange","purple","cyan4"))

#A few observations:
##Response variable histograms mostly vaguely normal (will be easy to standardize later)
##Pretty strong correlations exist among variables
##Likely easy to find principle components

#PCA analysis--------------------
peng.pca <- prcomp (~ bill_length_mm + bill_depth_mm + flipper_length_mm + body_mass_g,
                    data=penguins1,
                    scale. = TRUE) #Make sure to scale variables!

#Get factor loadings on principle components
peng.pca

#Visualize how much variation is explained by each principle component (Scree plot)
#Basic plot
plot(peng.pca)#however, this gives us the absolute variances, not % of total variance
#Use function from a different package (factoextra) for nicer plot
fviz_eig(peng.pca,addlabels = TRUE)
#PC1 and PC2 combined explain >88% of total variation
#Feel pretty good about plotting on the first 2 PCs

#Visualize results in a biplot--------
fviz_pca_biplot(peng.pca, label = "var", 
                col.ind = penguins1$species, palette = c("darkorange","purple","cyan4"), 
                col.var = "black", repel = TRUE,
                legend.title = "Species")
