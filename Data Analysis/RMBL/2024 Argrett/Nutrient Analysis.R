#Belowground Analysis
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2024 Argrett")

library(tidyverse)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plotting
library(car)#for regression analysis
library(emmeans)#post-hoc analysis
library(lme4)#for modeling linear mixed effect models
library(GGally)#this is new
library(factoextra)#this is new

soil <- read.csv("Soil Nutrients - Full.csv")
soil <- as.data.frame(unclass(soil),stringsAsFactors=TRUE)

soil <- plyr::rename(soil, c("nitrate" = "NO3",
                             "ammonium" = "NH4"))



nitrate <- ggplot(soil, aes(x = litter, y = K)) +
  geom_point(aes(color = (removal))) +
  facet_wrap(~burial) +
  labs(x = "Litter Treatment", y = "K")

nitrate

K.lmm <- lmer(K ~ litter*removal + burial + (1|block) + (1|pair), data = soil)
summary(K.lmm)
Anova(K.lmm)
emmip(K.lmm, litter ~ removal)
emmeans(K.lmm, pairwise ~  removal|litter)


#remove rows if NAs occur in any of the variables we care about
winter.soil <- subset(soil, burial == 'overwinter')
winter.soil <- subset( winter.soil, select = -c(B, S, Pb, Al, Cd))

winter.soil <- winter.soil %>% filter_at(vars(NO3, NH4, Ca, Mg, K, P, Fe, Mn, Cu, Zn),
                                    all_vars(!is.na(.)))

#Examine correlations among the body shape response variables------
#ggpairs plots is a good way to get lots of info (GGally package)
ggpairs(winter.soil, mapping = aes(color = litter), 
        columns = c("NO3", "NH4", 
                    "K", "P"))+
  scale_colour_manual(values = c("coral4","darkolivegreen3","deepskyblue4","darkgoldenrod2")) +
  scale_fill_manual(values = c("coral4","darkolivegreen3","deepskyblue4","darkgoldenrod2"))

#A few observations:
##Response variable histograms mostly vaguely normal (will be easy to standardize later)
##Pretty strong correlations exist among variables
##Likely easy to find principle components

#PCA analysis--------------------
winter.pca <- prcomp (~ NO3 + NH4 + Ca + Mg + K + P + Fe + Mn + Cu + Zn,
                    data=winter.soil,
                    scale. = TRUE) #Make sure to scale variables!

#Get factor loadings on principle components
winter.pca

#Visualize how much variation is explained by each principle component (Scree plot)
#Basic plot
plot(winter.pca)#however, this gives us the absolute variances, not % of total variance
#Use function from a different package (factoextra) for nicer plot
fviz_eig(winter.pca,addlabels = TRUE)
#PC1 and PC2 combined explain >88% of total variation
#Feel pretty good about plotting on the first 2 PCs

#Visualize results in a biplot--------
fviz_pca_biplot(winter.pca, label = "var", 
                col.ind = winter.soil$litter, palette = c("coral4","darkolivegreen3","deepskyblue4","darkgoldenrod2"), 
                col.var = "black", repel = TRUE,
                legend.title = "Litter")
