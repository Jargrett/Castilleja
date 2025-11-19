setwd("~/Desktop/Castilleja/Data Analysis/RMBL/Emerald Lake/Plant Data")
#---------------Data importing, cleaning, and resctructuring---------------#
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(conflicted)#helps reslove errors for similar functions between packages
library(rstatix)

#Declare conflicts
conflicts_prefer(rstatix::filter)

#load in cover data
cover.pre <- read.csv("Emerald Lake Plant Data - pre.csv")
cover.23 <- read.csv("Emerald Lake Plant Data - 2023.csv")
cover.24 <- read.csv("Emerald Lake Plant Data - 2024.csv")
cover.25 <- read.csv("Emerald Lake Plant Data - 2025.csv")
#combine datasets
cover.comb <- rbind.fill(cover.pre,cover.23,cover.24,cover.25)
cover.comb <- as.data.frame(unclass(cover.comb),stringsAsFactors=TRUE)
#remove castilleja and environmental rows for analysis
cover.comb.clean <- cover.comb[!(cover.comb$functional_group %in% "environmental"),]
cover.comb.clean <- cover.comb.clean[!(cover.comb.clean$code %in% "CASE"),]
comb.cov <- subset(cover.comb.clean, select = c('year','plot','code','percent_cover'))
#filter for year/pre and calculate
emerald.pre <- comb.cov %>% 
  filter (year == "Pre") %>%
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

m = list(emerald.23.matrix, emerald.24.matrix, emerald.25.matrix)
emerald.matrix <- bind_rows(m)
emerald.matrix %<>%  replace(is.na(.), 0)


#-------------------------Multivariate analysis-------------------------#
library(ggrepel)
library(vegan)
library(ggordiplots)
#plot data
setwd("~/Desktop/Castilleja/Data Analysis/RMBL/Emerald Lake")
plot <- read.csv("Emerald Lake Plot Data - Info.csv") #importing metadata
setwd("~/Desktop/Castilleja/Data Analysis/RMBL/Emerald Lake/Plant Data")
#calculate distance matrix
dist.23 <-vegdist(emerald.23.matrix, method="bray")
beta.23 <- betadisper(dist.23, plot$removal)
dist.25 <-vegdist(emerald.25.matrix, method="bray")
beta.25 <- betadisper(dist.25, plot$removal)
set.seed(20)

#Run NMDS on distance matrix
nmds.23 <- metaMDS(dist.23, distance="bray", #use bray-curtis distance
                k=3, #2 dimensions
                try=500) #for publication I recommend 500)
nmds.23#stress value 0.14 which is below .2 so we need to investigate


ordiplot(nmds.23, type="text", display="sites")

nmds.scores <- as.data.frame(vegan::scores(nmds.23))

NMDS <- cbind(plot,nmds.scores) #final dataset

perm <- adonis2(dist ~ removal*litter, data = NMDS, permutations=9999)
perm

ggplot(NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=removal , shape = litter), size = 2.2, alpha = 0.8) +
  scale_color_manual(values=c("#dda15e", "#606c38")) +
  coord_equal() +
  theme_bw()

