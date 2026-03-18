setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL Castilleja Observational Project")

#--------------------------Multivariate analysis-------------------------------#
#we will now run a (Multivariate analysis)
#This allows us to look at the compositional differences between our sites,castilleja,etc.
#we are working towards Matrix format so we can take our castilleja matrix as our starting point
library(ggrepel)
library(vegan)
library(ggordiplots)
library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(lme4)#for modeling linear mixed effect models
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(performance)#this is new
library(see)#this is new
library(lmerTest)
library(patchwork)
library(ggpubr)
library(rstatix)
library(permute)
library(sjPlot)
library(ggpmisc)

castilleja.cover <- readRDS("Processed Data/Total Castilleja Cover.rds")
species.matrix <- castilleja.cover[ -c(1:16)]
write.csv(species.matrix,"Processed Data/Species Matrix.csv", row.names = FALSE)
species.env <- subset(castilleja.cover, select=c(1:8))
write.csv(species.env,"Processed Data/Environmental Matrix.csv", row.names = FALSE)
#First calculate distance matrix
dist <-vegdist(species.matrix, method = "bray")

set.seed(20)
#Run NMDS on distance matrix
nmds <- metaMDS(dist, distance="bray", #use bray-curtis distance
                k=2, #2 dimensions
                try=500) #for publication I recommend 500)
nmds#stress value 0.14 which is below .2 so we need to investigate


ordiplot(nmds, type="text", display="sites")

nmds.scores <- as.data.frame(vegan::scores(nmds))

NMDS <- cbind(species.env,nmds.scores) #final dataset
saveRDS(NMDS,"NMDS.rds")

perm <- adonis2(dist ~ castilleja*species + castilleja*year + castilleja*site, 
                data = NMDS, permutations = 9999, by = "terms")
perm


cap.mod <- capscale( dist ~ castilleja*species + castilleja*year + castilleja*site + Condition(pair), 
                    data = NMDS)

perms <- how(blocks = NMDS$pair, nperm = 9999)

anova_res <- anova(cap.mod, permutations = perms, by = "terms")

capscale_table <- function(model, perms){
  
  a <- anova(model, permutations = perms, by = "terms")
  
  df <- data.frame(term = rownames(a), as.data.frame(a), 
                   row.names = NULL, check.names = FALSE)
  
  total_SS <- sum(df$SumOfSqs, na.rm = TRUE)
  df$R2 <- df$SumOfSqs / total_SS
  names(df)[names(df) == "Df"] <- "df"
  names(df)[names(df) == "SumOfSqs"] <- "SS"
  names(df)[names(df) == "F"] <- "F_value"
  names(df)[names(df) == "Pr(>F)"] <- "p_value"
  
  df
}

cap.table <- capscale_table(cap.mod, perms)
write.csv(cap.table, "Processed Data/Capscale Summary Output.csv", row.names=FALSE)

#
sites <- unique(NMDS$site)

results <- lapply(sites, function(s) {
  sub_data <- NMDS[NMDS$site == s, ]
  sub_dist <- as.dist(as.matrix(dist)[NMDS$site == s, NMDS$site == s])
  perms_sub <- how(blocks = sub_data$pair, nperm = 9999)
  cap_site <- capscale(sub_dist ~ castilleja + Condition(pair), data = sub_data)
  anov <- anova(cap_site, permutations = perms_sub)
  data.frame(site    = s, F_value = anov$F[1], p_value = anov$`Pr(>F)`[1])})

site_posthoc <- do.call(rbind, results)
site_posthoc$p_adj_bon <- p.adjust(site_posthoc$p_value, method = "bonferroni")
write.csv(site_posthoc, "Processed Data/Capscale Posthoc.csv", row.names=FALSE)


