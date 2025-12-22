setwd("~/Desktop/Castilleja/Data Analysis/RMBL Castilleja Experimental Project")
#----------Data importing, cleaning, and resctructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(conflicted)#helps reslove errors for similar functions between packages
library(car)

#Specifying conflicts
conflicted::conflicts_prefer(dplyr::recode)
conflicts_prefer(plyr::mutate)
conflicts_prefer(dplyr::filter)

#important cover data (raw)
cover.pre <- read.csv("Raw Data/Emerald Lake Plant Data - pre.csv")
cover.23 <- read.csv("Raw Data/Emerald Lake Plant Data - 2023.csv")
cover.24 <- read.csv("Raw Data/Emerald Lake Plant Data - 2024.csv")
cover.25 <- read.csv("Raw Data/Emerald Lake Plant Data - 2025.csv")

#combine datasets
cover.comb <- rbind.fill(cover.pre,cover.23,cover.24,cover.25)
cover.comb <- as.data.frame(unclass(cover.comb),stringsAsFactors=TRUE)
cover.comb %<>%
  mutate(removal = recode(removal,
                          "C" = "Control",
                          "R" = "Removal"))
#remove castilleja and environmental rows for analysis
cover.comb.clean <- cover.comb[!(cover.comb$functional_group %in% "environmental"),]
cover.comb.clean <- cover.comb.clean[!(cover.comb.clean$code %in% "CASE"),]
comb.cov <- subset(cover.comb.clean, select = c('year','plot','code','percent_cover'))
#filter for year/pre and calculate
emerald.pre <- comb.cov %>% 
  filter(year == "Pre") %>%
  select(-c(year))
emerald.23 <- comb.cov %>% 
  filter (year == "2023") %>%
  select(-c(year))
emerald.24 <- comb.cov %>% 
  filter (year == "2024") %>%
  select(-c(year))
emerald.25 <- comb.cov %>% 
  filter (year == "2025") %>%
  select(-c(year))

#convert to matrix format for diversity calculations
library(labdsv)#enables restructuring for ecological analysis
emerald.pre.matrix <- matrify(emerald.pre)
emerald.23.matrix <- matrify(emerald.23)
emerald.24.matrix <- matrify(emerald.24)
emerald.25.matrix <- matrify(emerald.25)

#---------------Diversity Calculations---------------#
library(vegan)#for calculating diversity
# Calculating Shannon diversity for plots
div.pre <- diversity(emerald.pre.matrix, index = "shannon")
div.23 <- diversity(emerald.23.matrix, index = "shannon")
div.24 <- diversity(emerald.24.matrix, index = "shannon")
div.25 <- diversity(emerald.25.matrix, index = "shannon")
# Calculating species richness for plots
rich.pre <- specnumber(emerald.pre.matrix)
rich.23 <- specnumber(emerald.23.matrix)
rich.24 <- specnumber(emerald.24.matrix)
rich.25 <- specnumber(emerald.25.matrix)
# Calculating species evenness for plots 
even.pre <- diversity(emerald.pre.matrix, index = "shannon") / log(specnumber(emerald.pre.matrix))
even.23 <- diversity(emerald.23.matrix, index = "shannon") / log(specnumber(emerald.23.matrix))
even.24 <- diversity(emerald.24.matrix, index = "shannon") / log(specnumber(emerald.24.matrix))
even.25 <- diversity(emerald.25.matrix, index = "shannon") / log(specnumber(emerald.25.matrix))

#---------------Combining results and exporting---------------#
plot <- read.csv("Raw Data/Emerald Lake Plot Data - Info.csv")
plot %<>%
  mutate(removal = recode(removal,
                          "C" = "Control",
                          "R" = "Removal"))

el.pre <- cbind(plot,div.pre,rich.pre,even.pre)
el.pre %<>% 
  plyr::mutate(year = 'pre') %>% 
  relocate(year) %>% 
  dplyr::rename(div = div.pre, rich = rich.pre, even = even.pre)
el.23 <- cbind(plot,div.23,rich.23,even.23)
el.23 %<>% 
  plyr::mutate(year = '2023') %>% 
  relocate(year) %>% 
  dplyr::rename(div = div.23, rich = rich.23, even = even.23)
el.24 <- cbind(plot,div.24,rich.24,even.24)
el.24 %<>% 
  plyr::mutate(year = '2024') %>% 
  relocate(year) %>% 
  dplyr::rename(div = div.24, rich = rich.24, even = even.24)
el.25 <- cbind(plot,div.25,rich.25,even.25)
el.25 %<>% 
  plyr::mutate(year = '2025') %>% 
  relocate(year) %>% 
  dplyr::rename(div = div.25, rich = rich.25, even = even.25)

diversity.pre <- rbind.fill(el.pre,el.23,el.24,el.25)
diversity <- rbind.fill(el.23,el.24,el.25)
write.csv(diversity, "Processed Data/Plant Diversity.csv", row.names=FALSE)

#----------Rank Abundance Curve Analysis----------#
library(codyn)
rac.cover <- cover.comb.clean %>% 
  filter (year != "Pre") %>% 
  subset(select = c('year','plot','pair','block','removal','litter','functional_group', 'code','count', 'percent_cover'))
rac.cover$year <- as.integer(rac.cover$year)

rac.diff <- RAC_difference(
  df = rac.cover,
  time.var = "year",
  species.var = "code",
  abundance.var = "percent_cover",
  replicate.var = "plot",
  treatment.var = "removal",
  block.var = "pair")

rac.change <- RAC_change(
  df = rac.cover,
  time.var = "year",
  species.var = "code",
  abundance.var = "percent_cover",
  replicate.var = "plot",
  reference.time = "1")

rac.change %<>% 
  filter (year2 != "2")

merged.rac.change <- left_join(plot, rac.change, by = "plot")

rs <- rank_shift(
  df = rac.cover,
  time.var = "year",
  species.var = "code",
  abundance.var = "percent_cover",
  replicate.var = "plot")


cs <- community_stability(rac.cover,
                    time.var = "year",
                    abundance.var = "percent_cover",
                    replicate.var = "plot")
turn <- turnover(
  df = rac.cover,
  time.var = "year",
  species.var = "code",
  abundance.var = "percent_cover",
  replicate.var = "plot",
  metric = "total")


#----------Analysis----------#
library(statmod)
library(lme4)
library(emmeans)#post-hoc analysis
library(car)#for regression analysis

rank.lm<- lmer(rank_change ~ litter*removal + (1|block), data = merged.rac.change)
summary(rank.lm)
Anova(rank.lm)
emmip(rank.lm, litter ~ removal)
emmeans(rank.lm, pairwise ~ removal|litter)
