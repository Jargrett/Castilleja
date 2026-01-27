setwd("~/Desktop/Castilleja/Data Analysis/RMBL Castilleja Experimental Project")
#----------Data importing, cleaning, and restructuring----------#
library(tidyverse)#for data wrangling and restructuring
library(magrittr)#for data wrangling and restructuring
library(plyr)#for data wrangling and restructuring
library(conflicted)#helps resolve errors for similar functions between packages
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

diversity.pre %<>% 
  mutate(year = recode(year,
                       "pre" = "0",
                       "2023" = "1",
                       "2024" = "2",
                       "2025" = "3",)) %>% 
  mutate(removal = recode(removal,
                          Control = "Present",
                          Removal = "Removed"))

#richness
rich.lmm <- lmer(rich ~ removal*year + (1|block) + (1|pair), data = diversity.pre)
summary(rich.lmm)
Anova(rich.lmm)#:removal:year p = 0.0008, Chisq = 16.8171, df = 3
emmeans(rich.lmm, pairwise ~ removal|year)
emmip(rich.lmm, removal ~ year)
library(ggeffects)

div.mean <- diversity.pre %>% 
  group_by(year,removal) %>% 
  dplyr::summarise(mean = mean(div),
                   se = sd(div)/sqrt(n()))

rich.mean <- diversity.pre %>% 
  group_by(year,removal) %>% 
  dplyr::summarise(mean = mean(rich),
                   se = sd(rich)/sqrt(n()))

even.mean <- diversity.pre %>% 
  group_by(year,removal) %>% 
  dplyr::summarise(mean = mean(even),
                   se = sd(even)/sqrt(n()))

div.plot <- ggplot(data = rich.mean, aes(x = year, y = mean, color = removal, group = removal)) +
  geom_point(shape=18, size = 4,position =  position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  geom_line(position =  position_dodge(width = 0.5)) +
  theme_pubr() +
  scale_color_manual(values = c("#333d29", "#b6ad90")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray", linewidth = 0.12)) +
  labs(x = "Growing Season", y = "Species Richness") 
div.plot

#-------Bare Ground Analysis------#
envi.cover <- cover.comb %>% 
  filter(functional_group == 'environmental') %>% 
  subset(select = c('year','plot','pair','block','removal','litter','code','percent_cover')) %>% 
  mutate(year = recode(year,
                       "Pre" = "0",
                       "2023" = "1",
                       "2024" = "2",
                       "2025" = "3",))
envi.cover$pair <- as.factor(envi.cover$pair)
envi.cover$plot <- as.factor(envi.cover$plot)
envi.cover$block <- as.factor(envi.cover$block)

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
split_and_name(envi.cover, "code")

total.envi <- envi.cover %>%
  group_by(year,plot,pair,block,removal,litter) %>%
  summarise(
    total_cover = sum(percent_cover, na.rm = TRUE)
  )

#bare ground
bare.lmm <- lmer(percent_cover ~ litter*removal + (1|block) + (1|pair), data = envi.cover_bare_X2025)
summary(bare.lmm)
Anova(bare.lmm)#(Higher in Removal): p = 0.0032, Chisq = 8.6926, df = 1
emmeans(bare.lmm, pairwise ~ litter|removal)
emmip(bare.lmm, litter ~ removal)

#litter ground
litter.lmm <- lmer(percent_cover ~ litter*removal + (1|block) + (1|pair), data = envi.cover_litter)
summary(litter.lmm)
Anova(litter.lmm)#No significant difference
emmeans(litter.lmm, pairwise ~ litter|removal)
emmip(litter.lmm, litter ~ removal)
#rock
rock.lmm <- lmer(percent_cover ~ litter*removal + (1|block) + (1|pair), data = envi.cover_rock)
summary(rock.lmm)
Anova(rock.lmm)#No significant difference
emmeans(rock.lmm, pairwise ~ litter|removal)
emmip(litter.lmm, litter ~ removal)

#Total Environment
total.lmm <- lmer(total_cover ~ year*removal + (1|year) + (1|block) + (1|pair), data = total.envi)
summary(total.lmm)
Anova(total.lmm)# removal = 
emmeans(total.lmm, pairwise ~ year|removal)
emmip(total.lmm, year ~ removal)


envi.total <- total.envi %>% 
  group_by(year,removal) %>% 
  dplyr::summarise(mean = mean(total_cover),
                   se = sd(total_cover)/sqrt(n()))

envi.total %<>% 
  dplyr::mutate(
    pat = ifelse(removal == "Removal", "stripe", "none")
  )

ggplot(envi.total, aes(x = removal, y = mean, fill = removal, pattern = pat)) +
  geom_bar_pattern(stat = "identity", color = "black", alpha = 0.8, width = 0.92, 
                   pattern_angle = 45, pattern_density = 0.12, 
                   pattern_spacing = 0.02, pattern_fill = '#333d29', pattern_colour = NA) +
  geom_errorbar(aes(ymin = mean, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("#333d29", "#b6ad90")) +
  scale_pattern_manual(values = c("none", "stripe")) +
  theme_pubr() +
  theme(legend.position = "none",panel.grid = element_blank()) +
  labs(x = "Parasite", y = "Environmental Cover (Bareground + Litter + Rock)")

#---------delta diversity Analysis---------#
library(statmod)
library(lme4)
library(emmeans) # for comparison of means
library(ggcharts)
library(ggthemes)
delta.div <- read.csv("Processed Data/Site Level Data - delta.csv")
conflicts_prefer(lme4::lmer)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::arrange)

div.lmm <- lmer(delta_rich ~ litter*removal + (1|block) + (1|pair), data = delta.div)
summary(div.lmm)
Anova(div.lmm)
emmip(div.lmm, litter~removal)
emmeans(div.lmm, pairwise ~ litter|removal)

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

arrow_info <- delta %>%
  arrange(plot,year) %>%
  group_by(plot) %>%
  mutate(rich_lag = lag(rich),
         rich_diff = rich - rich_lag) %>%
  select(plot,rich_diff) %>%
  filter(!is.na(rich_diff)) %>%
  mutate(arrow_draw = ifelse(rich_diff==0,"no arrow","arrow"))

delta <- merge(delta,arrow_info,by="plot")

delta_plot <- delta %>%
  arrange(plot, year) %>%
  group_by(plot) %>%
  mutate(
    field_plot2 = as.numeric(gsub("A|B", "", field_plot)),
    rich_max = max(rich),
    rich_2023 = lag(rich)
  ) %>%
  ungroup() %>%
  mutate(
    text_nudge_x = case_when(rich == rich_max ~ 1.2,
                             rich <= rich_max ~ -1.2,
                             rich >= rich_max ~ 1.2)
  )

delta_plot %>% 
  ggplot(aes(x = rich,y = field_plot2)) +
  geom_path(aes(group = plot), color = "#b7b7a4", linewidth = 0.5) +
  geom_segment(
    data = . %>% filter(year == 2025 & arrow_draw == "arrow"),
    aes(x = rich_2023, xend = ifelse(rich_diff > 0, rich - 0.55, rich + 0.55)),
    arrow = arrow(angle = 25, type = "closed", length = unit(0.24, "cm")),
    color = "gray55") +
  geom_point(aes(color=as.factor(year)), size=3.5) +
  geom_text(aes(label = rich), size = 3.25, nudge_x = delta_plot$text_nudge_x) +
  scale_color_manual(values = c("#b6ad90", "#656d4a")) +
  theme_classic() +
  labs_pubr() + 
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray23", linewidth = 0.12)) +
  xlim(6,26) +
  scale_y_continuous(breaks=seq(1,20,by=1))+
  facet_wrap(~removal)+
  labs(x = "Species Richness", y = "Paired Plot")

#----------Nearest Neighbor Analysis----------#
comb.cov.NN <- subset(cover.comb.clean, select = c('year','plot','code','functional_group',
                                             'percent_cover','nearest_neighbor'))
write.csv(comb.cov.NN, "Processed Data/Neighbor Cover.csv", row.names=FALSE)
  

