setwd("~/Desktop/Castilleja/Data Analysis/RMBL Castilleja Experimental Project")
#----------Data importing, cleaning, and resctructuring----------#
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
library(magrittr)#for data wrangling and restructuring

soil <- read.csv("Raw Data/Soil Nutrients - Full.csv")
soil.ex <- read.csv("Raw Data/Soil Nutrients - Full Excluded.csv")
soil <- as.data.frame(unclass(soil),stringsAsFactors=TRUE)

soil.full <- filter(soil, burial == "full study")
soil.overwinter <- filter(soil, burial == "overwinter")
soil.within <- filter(soil, burial == "within")
soil.year <- filter(soil, burial == "year")

#-------PCA Analysis-------#
library(stats)
library(vegan)
library(tidyverse)
library(ggbiplot)
library(GGally)#this is new
library(factoextra)#this is new
# Units: (micro grams/10cm2/burial length)


#Short version
#PCA analysis--------------------

full.pca <- prcomp (~ NO3 + NH4 + Ca + Mg + K + P + Fe + Mn + Cu + Zn + B + S + Pb + Al + Cd,
                    data=soil.full,
                    scale. = TRUE)

#Get factor loadings on principle components
full.pca

#Visualize how much variation is explained by each principle component (Scree plot)
#Basic plot
plot(full.pca)#however, this gives us the absolute variances, not % of total variance
#Use function from a different package (factoextra) for nicer plot
fviz_eig(full.pca,addlabels = TRUE)

pc_scores <- full.pca$x  #
#PC1 and PC2 combined explain >55% of total variation
#Feel pretty good about plotting on the first 2 PCs
dist_matrix <- dist(pc_scores[, 1:3], method = "euclidean")

full.perm <- adonis2(dist_matrix ~ litter*removal, data = soil.full, permutations = 9999)
print(full.perm)

#Visualize results in a biplot--------
fviz_pca_biplot(full.pca, label = "var",
                addEllipses = TRUE,
                ellipse.level=0.95,
                habillage=soil.full$litter,
                col.ind = full.pca$litter, palette = c("#4B3B40","#9B7E46","#B3B1BE","#909256"), 
                col.var = "black", repel = TRUE,
                legend.title = "Litter Treatment") +
  theme_minimal()


#------ NPKCa Analysis-------#
#Nitrate
NO3.lme <- lmer(NO3 ~ litter*removal + (1|block) + (1|pair), data = soil.full)
summary(NO3.lme)
Anova(NO3.lme) #No significance
emmip(NO3.lme, litter ~ removal)
emmeans(NO3.lme, pairwise ~  removal|litter)

#Ammonium
NH4.lme <- lmer(NH4 ~ litter*removal + (1|block) + (1|pair), data = soil.full)
summary(NH4.lme)
Anova(NH4.lme) #No significance, litter p = 0.087
emmip(NH4.lme, litter ~ removal)
emmeans(NH4.lme, pairwise ~  removal|litter)

#Phosphorus
P.lme <- lmer(P ~ litter*removal + (1|block) + (1|pair), data = soil.full)
summary(P.lme)
Anova(P.lme) #No significance
emmip(P.lme, litter ~ removal)
emmeans(P.lme, pairwise ~  removal|litter)

#Potassium
K.lme <- lmer(K ~ litter*removal + (1|block) + (1|pair), data = soil.full)
summary(K.lme)
Anova(K.lme) #litter: Chisq = 20.789, p = 0.0001 *, removal: Chisq = 3.1979, p = 0.074.
emmip(K.lme, litter ~ removal)
emmeans(K.lme, pairwise ~  removal|litter)

#Calcium
Ca.lme <- lmer(Ca ~ litter*removal + (1|block) + (1|pair), data = soil.full)
summary(Ca.lme)
Anova(Ca.lme) #No significance
emmip(Ca.lme, litter ~ removal)
emmeans(Ca.lme, pairwise ~  removal|litter)



#Nutrient graphs
soil.long <- soil.full %>%
  select(-17:-25) %>% 
  pivot_longer(
    cols = 11:16 ,      # columns to pivot
    names_to = "nutrient",  # new column for former col names
    values_to = "value"        # new column for the numeric values
  )

soil.mean <- soil.long %>% 
  group_by(nutrient,litter,removal) %>% 
  dplyr::summarise(mean = mean(value),
                   se = sd(value)/sqrt(n()))

soil.mean %<>% 
  dplyr::mutate(
    pat = ifelse(removal == "R", "stripe", "none"))

split_and_name <- function(df, column) {
  # Get the dataframe name as a string
  df_name <- deparse(substitute(df))
  # Ensure column exists
  if (!column %in% colnames(df)) {
    stop("Column not found in dataframe")
  }
  # Split the dataframe by the column
  df_list <- split(df, df[[column]])
  # Remove dataframes with 0 rows
  df_list <- df_list[sapply(df_list, nrow) > 0]
  # If nothing remains, stop gracefully
  if (length(df_list) == 0) {
    warning("All split dataframes have 0 rows; nothing to assign.")
    return(invisible(NULL))
  }
  # Make valid R names for safety
  names(df_list) <- paste0(df_name, "_", make.names(names(df_list)))
  # Assign each dataframe into the global environment
  list2env(df_list, envir = .GlobalEnv)
  # Return the list invisibly
  invisible(df_list)
}

split_and_name(soil.mean, "nutrient") 

nitrate <- ggplot(soil.mean_NO3, aes(x = litter, y = mean, fill = removal, pattern = pat)) +
  geom_bar_pattern(stat = "identity", color = "black", alpha = 0.8, width = 0.92,
                   position = position_dodge(width = 0.92),
                   pattern_angle = 45, pattern_density = 0.12, 
                   pattern_spacing = 0.02, pattern_fill = '#333d29', pattern_colour = NA) +
  geom_errorbar(aes(ymin = mean, ymax = mean + se), width = 0.2, position = position_dodge(width = 0.92)) +
  scale_fill_manual(values = c("#333d29","#b6ad90")) +
  scale_pattern_manual(values = c("none", "stripe")) +
  theme_pubr() +
  theme(legend.position = "none",panel.grid = element_blank()) +
  labs(x = "Litter Treatment", y = "NO3 (mg/10cm2/t)")
nitrate

ammonium <- ggplot(soil.mean_NH4, aes(x = litter, y = mean, fill = removal, pattern = pat)) +
  geom_bar_pattern(stat = "identity", color = "black", alpha = 0.8, width = 0.92,
                   position = position_dodge(width = 0.92),
                   pattern_angle = 45, pattern_density = 0.12, 
                   pattern_spacing = 0.02, pattern_fill = '#333d29', pattern_colour = NA) +
  geom_errorbar(aes(ymin = mean, ymax = mean + se), width = 0.2, position = position_dodge(width = 0.92)) +
  scale_fill_manual(values = c("#333d29","#b6ad90")) +
  scale_pattern_manual(values = c("none", "stripe")) +
  theme_pubr() +
  theme(legend.position = "none",panel.grid = element_blank()) +
  labs(x = "Litter Treatment", y = "NH4 (mg/10cm2/t)")
ammonium

phosphorus <- ggplot(soil.mean_P, aes(x = litter, y = mean, fill = removal, pattern = pat)) +
  geom_bar_pattern(stat = "identity", color = "black", alpha = 0.8, width = 0.92,
                   position = position_dodge(width = 0.92),
                   pattern_angle = 45, pattern_density = 0.12, 
                   pattern_spacing = 0.02, pattern_fill = '#333d29', pattern_colour = NA) +
  geom_errorbar(aes(ymin = mean, ymax = mean + se), width = 0.2, position = position_dodge(width = 0.92)) +
  scale_fill_manual(values = c("#333d29","#b6ad90")) +
  scale_pattern_manual(values = c("none", "stripe")) +
  theme_pubr() +
  theme(legend.position = "none",panel.grid = element_blank()) +
  labs(x = "Litter Treatment", y = "P (mg/10cm2/t)")
phosphorus

potassium <- ggplot(soil.mean_K, aes(x = litter, y = mean, fill = removal, pattern = pat)) +
  geom_bar_pattern(stat = "identity", color = "black", alpha = 0.8, width = 0.92,
                   position = position_dodge(width = 0.92),
                   pattern_angle = 45, pattern_density = 0.12, 
                   pattern_spacing = 0.02, pattern_fill = '#333d29', pattern_colour = NA) +
  geom_errorbar(aes(ymin = mean, ymax = mean + se), width = 0.2, position = position_dodge(width = 0.92)) +
  scale_fill_manual(values = c("#333d29","#b6ad90")) +
  scale_pattern_manual(values = c("none", "stripe")) +
  theme_pubr() +
  theme(legend.position = "none",panel.grid = element_blank()) +
  labs(x = "Litter Treatment", y = "K (mg/10cm2/t)")
potassium

nutrientplots <- ggarrange(nitrate, ammonium, phosphorus, potassium,
                          labels = c("A", "B","C", "D"), 
                          nrow = 2, ncol = 2)

nutrientplots 
