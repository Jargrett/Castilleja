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



NMDS <- readRDS("Processed Data/Community NMDS.rds")
dist.all <- readRDS("Processed Data/dist_all_bray.rds")

bd_all <- betadisper(dist.all, NMDS$removal)
permutest(bd_all, permutations = 999)
boxplot(bd_all)

cap.mod <- capscale(dist.all ~ removal * litter * year + Condition(pair),
                    data = NMDS)

perms <- how(blocks = plot_meta$pair, nperm = 9999)
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
