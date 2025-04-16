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
#PC1 and PC2 combined explain >55% of total variation
#Feel pretty good about plotting on the first 2 PCs

#Visualize results in a biplot--------
fviz_pca_biplot(winter.pca, label = "var",
                col.ind = winter.soil$litter, palette = c("coral4","darkgoldenrod2","deepskyblue4","darkolivegreen3"), 
                col.var = "black", repel = TRUE,
                legend.title = "Litter Treatment")

perm <- adonis2(winter.pca ~ litter*removal, data = NMDS, permutations=9999)
perm


#-------PCA Analysis-------#
library(stats)
library(vegan)
library(tidyverse)
library(ggbiplot)
library(GGally)#this is new
library(factoextra)#this is new
# Units: (micro grams/10cm2/burial length)

#Create an object that contains the response variables (Nutrients)
# Units: (micro grams/10cm2/burial length)
nutrients <- subset( winter.soil, select = c(NO3, NH4, Ca, Mg, K, P, Fe, Mn, Cu, Zn))

#subset out predictors
envi <- subset( winter.soil, select = c(5:10))

#Normalize data – subtract the mean and divide by the standard deviation for each column.
nutrients.scale <- scale(nutrients)

#Convert raw data to correlation matrix
nutrients.matrix <- cov(nutrients.scale)
nutrients.matrix %>% round(2)

#Conduct eigenanalysis
nutrients.eigen <- eigen(nutrients.matrix)

#There will be 10 eigenvalues (equal to number of response variables)
nutrients.eigen$values %>% round(3)

#The eigenvalues sum to the total of the diagonal of the matrix analyzed
sum(nutrients.eigen$values)

#Express eigenvalues as a proportion of the total amount of variation
nutrients.eigen.prop <- nutrients.eigen$values / sum(nutrients.eigen$values)
nutrients.eigen.prop %>% round(3)

#visualize eigenvectors
#rows correspond to the original variables
#columns correspond to the principal component identified by each eigenvalue.
nutrients.eigen$vectors %>% round(2)

#loadings are the coefficients of the linear equation that ‘blend’ the variables together for that principal component.
loadings <- nutrients.eigen$vectors[ , 1:3] %>%
data.frame(row.names = colnames(nutrients)) %>%
  dplyr::rename("PC1" = X1, "PC2" = X2, "PC3" = X3) %>%
  round(digits = 3)
loadings

#Interpret the Scores of Sample Units
pc.scores <- nutrients.scale %*% nutrients.eigen$vectors
pc.scores %>% head()

#The PC scores are orthogonal to (i.e., uncorrelated with) each other
cor(pc.scores) %>% round(3) 

nutrients.pca <- princomp(nutrients.scale, cor = TRUE)

scores <- nutrients.pca$scores

nutrients.bind <- cbind(envi,scores) 

summary(nutrients.pca, loadings = TRUE, cutoff = 0)

fviz_eig(nutrients.pca,addlabels = TRUE)


#-------#
adonis2(nutrients.bind$Comp.1 ~ litter*removal,
        data = nutrients.bind,
        method = "euc",
        permutations=9999)

ggplot(data = nutrients.bind, aes(x = litter, y = Comp.1)) +
  geom_boxplot() +
  geom_jitter(aes(colour = litter), width = 0.3, height = 0) +
  theme_bw()

