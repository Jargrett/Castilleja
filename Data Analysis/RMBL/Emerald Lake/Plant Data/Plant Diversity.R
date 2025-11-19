setwd("~/Desktop/Castilleja/Data Analysis/RMBL/Emerald Lake/Plant Data")

#---------------Data importing, cleaning, and resctructuring---------------#
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(conflicted)#helps reslove errors for similar functions between packages
library(car)

conflicted::conflicts_prefer(dplyr::recode)
conflicts_prefer(plyr::mutate)
conflicts_prefer(dplyr::filter)
#load in cover data
cover.pre <- read.csv("Emerald Lake Plant Data - pre.csv")
cover.23 <- read.csv("Emerald Lake Plant Data - 2023.csv")
cover.24 <- read.csv("Emerald Lake Plant Data - 2024.csv")
cover.25 <- read.csv("Emerald Lake Plant Data - 2025.csv")
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

m = list(emerald.23.matrix, el.24.matrix, el.25.matrix)
emerald.matrix <- bind_rows(m)
emerald.matrix %<>%  replace(is.na(.), 0) %>% 
  relocate(year)

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
setwd("~/Desktop/Castilleja/Data Analysis/RMBL/Emerald Lake")
plot <- read.csv("Emerald Lake Plot Data - Info.csv")
plot %<>%
  mutate(removal = recode(removal,
                          "C" = "Control",
                          "R" = "Removal"))
setwd("~/Desktop/Castilleja/Data Analysis/RMBL/Emerald Lake/Plant Data")
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
write.csv(diversity, "~/Desktop/Castilleja/Data Analysis/RMBL/Emerald Lake/Plant Data/Plant Diversity.csv", row.names=FALSE)

#---------------Delta Analysis and visualization---------------#
library(statmod)
library(lme4)
library(emmeans) # for comparison of means
library(ggcharts)
library(ggthemes)
delta.div <- read.csv("Site Level Data - delta.csv")
conflicts_prefer(lme4::lmer)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::arrange)

div.lmm <- lmer(delta_rich ~ litter*removal + (1|block) + (1|pair), data = delta.div)
summary(div.lmm)
Anova(div.lmm)
emmip(div.lmm, litter~removal)
emmeans(div.lmm, pairwise ~ removal|litter)

max.div <- diversity %>% 
  filter(year != "2024")%>% 
  group_by(plot,pair,removal) %>% 
  reframe(max = max(rich))
delta <- merge(max.div,diversity, by =c('plot'))

delta %<>% 
  filter(year != "2024")%>% 
  select(-c(pair.y,removal.y)) %>% 
  dplyr::rename(pair = pair.x,
                removal = removal.x)
nudge_value=.6
d.m <- delta %>% 
  mutate(field_plot2 = as.numeric(gsub("A|B","",field_plot))) %>% 
ggplot(aes(x = rich,y = field_plot2)) +
  geom_path(aes(group = plot), color="#b7b7a4", linewidth=0.5, arrow = 
              j(angle=10, type="closed",length=unit(0.5,"cm"))) + 
  geom_point(aes(color=year), size=3) +
  scale_color_manual(values=c("#dda15e", "#606c38")) +
  #theme_pubr() +
  theme_par()+
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  geom_text(aes(label=rich), size = 3.25,
            nudge_x=if_else(
              delta$rich==delta$max, nudge_value, -nudge_value), 
            hjust=if_else(delta$rich==delta$max,0,1)) +

  xlim(7,26) +
  scale_y_continuous(breaks=seq(1,20,by=1))+
  facet_wrap(~removal)+
  labs(x = "Species Richness", y = "Paired Plot")
d.m

#----------------Functional diversity Analysis----------------#
conflicts_prefer(dplyr::summarize)
conflicts_prefer(dplyr::c_across)
cover_sum <- cover.comb %>% 
  group_by(year,plot,functional_group) %>%
  filter (year != "Pre") %>%
  summarize(plot_cover = sum(percent_cover),
            removal = unique(removal),
            litter = unique(litter),
            block = unique(block),
            pair = unique(pair))

cover_sum_wide <- cover_sum %>% 
  pivot_wider(names_from = functional_group,
              values_from = plot_cover) %>% 
  replace(is.na(.), 0)  %>%
  mutate(total_cover = sum(c_across(environmental:shrub))) %>% 
  mutate(plant_cover = sum(c(forb, grass, legume, sedge, shrub)))

#Functional richness
func.cov <- subset(cover.comb.clean, select = c('year','plot','functional_group'))

func.pre <- func.cov %>% 
  filter (year == "Pre") %>%
  count(plot, functional_group, sort = TRUE)
func.23 <- func.cov %>% 
  filter (year == "2023") %>%
  count(plot, functional_group, sort = TRUE)
func.24 <- func.cov %>% 
  filter (year == "2024") %>%
  count(plot, functional_group, sort = TRUE)
func.25 <- func.cov %>% 
  filter (year == "2025") %>%
  count(plot, functional_group, sort = TRUE)

func.total <- func.cov %>% 
  count(year, plot, functional_group)


func.total.wide <- func.total %>% 
  pivot_wider(names_from = functional_group,
              values_from = n) %>% 
  replace(is.na(.), 0) 

write.csv(func.total.wide, "~/Desktop/Castilleja/Data Analysis/RMBL/Emerald Lake/Plant Data/Functional Diversity.csv", row.names=FALSE)

#
func.lmm <- lmer(delta_sedge_rich ~ removal*litter + (1|block) + (1|pair), data = delta.div)
summary(func.lmm)
Anova(func.lmm)
emmip(func.lmm, litter~removal)
emmeans(func.lmm, pairwise ~ removal|litter)
