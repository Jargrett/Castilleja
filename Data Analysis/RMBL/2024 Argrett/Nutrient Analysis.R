#Belowground Analysis
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2024 Argrett")

library(tidyverse)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(ggplot2)#for plotting
library(ggpubr)#extended functions for plottinglibrary(remotes)
library(ggpattern)
library(car)#for regression analysis
library(emmeans)#post-hoc analysis
library(lme4)#for modeling linear mixed effect models
library(GGally)#this is new
library(factoextra)#this is new

soil <- read.csv("Soil Nutrients - Full.csv")
soil.ex <- read.csv("Soil Nutrients - Full Excluded.csv")
soil <- as.data.frame(unclass(soil),stringsAsFactors=TRUE)

soil.full <- filter(soil, burial == "full study")
soil.overwinter <- filter(soil, burial == "overwinter")
soil.within <- filter(soil, burial == "within")
soil.year <- filter(soil, burial == "year")

nitrate <- ggplot(soil, aes(x = litter, y = NH4)) +
  geom_point(aes(color = (removal))) +
  facet_wrap(~burial) +
  labs(x = "Litter Treatment", y = "NO3")

nitrate


nitrate.lme <- lmer(K ~ litter*removal + (1|block) + (1|pair), data = soil.full)
summary(nitrate.lme)
Anova(nitrate.lme)
emmip(nitrate.lme, litter ~ removal)
emmeans(nitrate.lme, pairwise ~  removal|litter)

soil.NO3 <- soil.full %>% 
  group_by(litter,removal,burial) %>% 
  dplyr::summarise(mean= mean(K),
                   se = sd(K)/sqrt(n()))

NO3.plot <- ggplot(data = soil.NO3, aes(x = litter, y = mean, fill = removal, pattern = removal)) +
  geom_bar(stat="identity", color="black",position=position_dodge2(preserve="single", padding=0), alpha=0.5) +
  geom_errorbar(aes(ymin=mean, ymax=mean+se), alpha=1, width=.2,position=position_dodge(.9)) +
  theme(panel.grid = element_blank(),legend.position = "bottom") +
  labs(x="Litter Type", y="Nitrate (micro grams/10cm2/burial length)", pattern="") + scale_fill_manual(values=cols)+
  geom_bar_pattern(stat="identity", position = "dodge", color = "black", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  scale_fill_manual( values=c("#004D40", "#FFC107"))
NO3.plot

a<-ggplot(COroots, aes(x = Invasion, y = mean, fill = grp, pattern=Treatment, alpha=alfa)) +
  geom_bar(stat="identity", color="black",position=position_dodge2(preserve="single", padding=0), alpha=0.5) +
  geom_errorbar(aes(ymin=mean, ymax=mean+se), alpha=1, width=.2,position=position_dodge(.9)) +
  theme(panel.grid = element_blank(),legend.position = "bottom") +
  labs(x="Soil Type", y="Root Mass (g)", pattern="") + scale_fill_manual(values=cols)+
  geom_bar_pattern(stat="identity", position = "dodge", color = "black", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c("stripe", "none")) +
  theme(axis.text.x = element_text(size=11, colour= "black", angle = 0),axis.text.y = element_text(size=12, colour="black"),strip.text = element_text(size = 10))+
  theme(axis.title.x = element_text(size=12, colour="black"), axis.title.y = element_text(size= 12, colour = "black"))+
  theme(legend.title = element_blank())+theme(axis.ticks.x = element_blank())+
  theme(panel.grid = element_blank(),legend.position = "bottom")+
  theme(legend.title = element_blank())+theme(legend.text = element_text(size = 10, colour = "black"))+
  theme(legend.background = element_blank())+ylim(0, 2.0) +
  facet_grid(~ Grown, scales = "free") +
  geom_signif(data=data.frame(Grown=c("Linaria", "Native")), aes(xmin =c(1,1), xmax = c(2,2), y_position = c(0.95, 1.95), annotations=c("*","*")), tip_length = 0, manual=T, inherit.aes = FALSE) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),fill = guide_legend(override.aes = list(pattern = "none"))) +
  guides(fill=FALSE)+ guides(size="none", alpha="none")


#remove rows if NAs occur in any of the variables we care about
soil.full <- subset(soil.full, select = -c(B, S, Pb, Al, Cd))

soil.full <- soil.full %>% filter_at(vars(NO3, NH4, Ca, Mg, K, P, Fe, Mn, Cu, Zn),
                                    all_vars(!is.na(.)))

#Examine correlations among the body shape response variables------
#ggpairs plots is a good way to get lots of info (GGally package)
ggpairs(soil.full, mapping = aes(color = litter), 
        columns = c("NO3", "NH4", 
                    "K", "P"))+
  scale_colour_manual(values = c("coral4","darkolivegreen3","deepskyblue4","darkgoldenrod2")) +
  scale_fill_manual(values = c("coral4","darkolivegreen3","deepskyblue4","darkgoldenrod2"))

#A few observations:
##Response variable histograms mostly vaguely normal (will be easy to standardize later)
##Pretty strong correlations exist among variables
##Likely easy to find principle components

#PCA analysis--------------------
full.pca <- prcomp (~ NO3 + NH4 + Ca + Mg + K + P + Fe + Mn + Cu + Zn,
                    data=soil.full,
                    scale. = TRUE) #Make sure to scale variables!

#Get factor loadings on principle components
full.pca

#Visualize how much variation is explained by each principle component (Scree plot)
#Basic plot
plot(full.pca)#however, this gives us the absolute variances, not % of total variance
#Use function from a different package (factoextra) for nicer plot
fviz_eig(full.pca,addlabels = TRUE)
#PC1 and PC2 combined explain >55% of total variation
#Feel pretty good about plotting on the first 2 PCs

#Visualize results in a biplot--------
fviz_pca_biplot(full.pca, label = "var",
                col.ind = full.pca$litter, palette = c("coral4","darkgoldenrod2","deepskyblue4","darkolivegreen3"), 
                col.var = "black", repel = TRUE,
                legend.title = "Litter Treatment")

perm <- adonis2(full.pca ~ litter*removal, data = NMDS, permutations=9999)
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
year.soil <- subset(soil, burial == 'year')
year.soil <- subset( year.soil, select = -c(B, S, Pb, Al, Cd))

nutrients <- subset(soil.full, select = c(NO3, NH4, Ca, Mg, K, P, Fe, Mn, Cu, Zn))

#subset out predictors
envi <- subset(soil.full, select = c(5:10))

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

