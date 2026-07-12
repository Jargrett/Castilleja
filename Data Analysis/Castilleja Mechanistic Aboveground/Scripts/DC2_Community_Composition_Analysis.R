setwd("~/Desktop/Castilleja/Data Analysis/Castilleja Mechanistic Aboveground")
#----------Data importing, cleaning, and restructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(statmod)
library(lme4)#for modeling linear mixed effect models
library(emmeans)#post-hoc analysis
library(car)#for regression analysis
library(performance)#this is new
library(see)#this is new
library(lmerTest)
library(rstatix)
library(magrittr)#for data wrangling and restructuring
library(conflicted)
library(vegan)


NMDS <- readRDS("Processed Data/Community NMDS.rds")
dist.all <- readRDS("Processed Data/dist_all_bray.rds")
 
#----------------------------------------------------------#
#--------------MULTIVARIATE DISPERSION (PERMDISP)----------#
#----------------------------------------------------------#
# Run BEFORE interpreting capscale/PERMANOVA results.
# A significant PERMANOVA can reflect a shift in group centroids (a true
# compositional difference) OR unequal within-group dispersion. PERMDISP
# distinguishes these: a NON-significant result here means the compositional
# difference reflects a centroid shift, which is what we want to claim.
 
set.seed(20)
 
bd_removal <- betadisper(dist.all, NMDS$removal)
bd_litter  <- betadisper(dist.all, NMDS$litter)
bd_year    <- betadisper(dist.all, NMDS$year)
 
pt_removal <- permutest(bd_removal, permutations = 9999)
pt_litter  <- permutest(bd_litter,  permutations = 9999)
pt_year    <- permutest(bd_year,    permutations = 9999)
 
pt_removal
pt_litter
pt_year
 
# tidy summary table for reporting
disp_tab <- tibble::tibble(
  factor = c("Removal", "Litter", "Year"),
  df     = c(pt_removal$tab$Df[1], pt_litter$tab$Df[1], pt_year$tab$Df[1]),
  F      = c(pt_removal$tab$F[1],  pt_litter$tab$F[1],  pt_year$tab$F[1]),
  p      = c(pt_removal$tab$`Pr(>F)`[1],
             pt_litter$tab$`Pr(>F)`[1],
             pt_year$tab$`Pr(>F)`[1]))
disp_tab
 
write.csv(disp_tab, "Processed Data/PERMDISP.csv", row.names = FALSE)
 
# mean distance-to-centroid per group (useful for the results sentence)
tapply(bd_removal$distances, NMDS$removal, mean)
tapply(bd_litter$distances,  NMDS$litter,  mean)
tapply(bd_year$distances,    NMDS$year,    mean)
 
#----------------------------------------------------------#
#------------COMPOSITION: DISTANCE-BASED RDA---------------#
#----------------------------------------------------------#
cap.mod <- capscale(dist.all ~ removal * litter * year + Condition(pair),
                    data = NMDS)
 
perms <- how(blocks = NMDS$pair, nperm = 9999)
anova_res <- anova(cap.mod, permutations = perms, by = "terms")
anova_res
 
litters <- unique(NMDS$litter)
 
results <- lapply(litters, function(l) {
  sub_data <- droplevels(NMDS[NMDS$litter == l, ])
  sub_dist <- as.dist(as.matrix(dist.all)[NMDS$litter == l, NMDS$litter == l])
  perms_sub <- how(blocks = sub_data$pair, nperm = 9999)
  cap_litter <- capscale(sub_dist ~ removal + Condition(pair), data = sub_data)
  anov <- anova(cap_litter, permutations = perms_sub)
  data.frame(litter  = l,
             F_value = anov$F[1],
             p_value = anov$`Pr(>F)`[1])
})
 
litter_posthoc <- do.call(rbind, results)
litter_posthoc$p_adj_bon <- p.adjust(litter_posthoc$p_value, method = "bonferroni")
write.csv(litter_posthoc, "Processed Data/Capscale Posthoc.csv", row.names = FALSE)
 
