# Emerald Lake plant cover
# Initiated: 8/9/24
# Completed: TBD

#####################
#                   # - Cover analysis of paired plots
#      Question     # - Total Castilleja cover
#                   # - NN analysis and richness
##################### - 

# Set working directory

setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2024 Argrett")

# Load-in packages
library(tidyverse) # for data working
library(vegan) # for diversity
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car)#for regression analysis
library(ggpubr)#extended functions for plotting
library(emmeans) # for comparison of means
library(rstatix) # for comparison of means
library(labdsv)
library(emmeans)#post-hoc analysis
library(indicspecies)
library(plyr)
library(statmod)
library(rstatix)
library(lme4)

EL.23 <- read.csv("Emerald Lake Plant Data - 2023.csv")
EL.24 <- read.csv("Emerald Lake Plant Data - 2024.csv")

#-----------------data structure-----------------#
#Combining Dataset
EL.comb <- rbind.fill(EL.23,EL.24)
str(EL.comb)
summary(EL.comb)

#changing to factors
EL.comb <- as.data.frame(unclass(EL.comb),stringsAsFactors=TRUE)
#Removing pre treatment analysis from dataset
emerald <- EL.comb%>% filter (collection == "Post")
str(emerald)
summary(emerald)

#------------Diversity Analysis------------#
# here we will be working with the cover column
# we will need to conververt data to a matrix format
emerald.24 <- EL.comb%>% filter (year == "2024")
emerald.24 <- emerald.24[!(emerald.24$functional_group %in% "environmental"),]
emerald.24 <- emerald.24[!(emerald.24$code %in% "CASE"),]
emerald.24.cov <- subset(emerald.24, select = c('plot','code','cover'))

emerald.23 <- EL.comb%>% filter (year == "2023")
emerald.23 <- emerald.23[!(emerald.23$functional_group %in% "environmental"),]
emerald.23 <- emerald.23[!(emerald.23$code %in% "CASE"),]
emerald.23.cov <- subset(emerald.23, select = c('plot','code','cover'))

emerald.24.matrix <- matrify(emerald.24.cov)
emerald.23.matrix <- matrify(emerald.23.cov)
# Calculating Shannon diversity for plots
div.24 <- diversity(emerald.24.matrix, index = "shannon")
div.23 <- diversity(emerald.23.matrix, index = "shannon")
# Calculating species richness for plots
rich.24 <- specnumber(emerald.24.matrix)
rich.23 <- specnumber(emerald.23.matrix)
# Calculating species evenness for plots 
even.24 <- diversity(emerald.24.matrix, index = "shannon") / log(specnumber(emerald.24.matrix)) 
even.23 <- diversity(emerald.23.matrix, index = "shannon") / log(specnumber(emerald.23.matrix)) 

#combined data set with Plot Data File and calculated values
p.data <- read.csv("Emerald Lake Plot Data - Info.csv") #importing metadata
el.div.24 <- cbind(p.data,div.24,rich.24,even.24) #final dataset
el.div.23 <- cbind(p.data,div.23,rich.23,even.23)

el.div.24 <- plyr::rename(el.div.24, c("div.24" = "div",
                                       "even.24" = "even",
                                       "rich.24" = "rich"))

el.div.23 <- plyr::rename(el.div.23, c("div.23" = "div",
                                       "even.23" = "even",
                                       "rich.23" = "rich"))

write.csv(el.div.24, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2024 Argrett/Diversity 2024.csv", row.names=FALSE)
write.csv(el.div.23, "/Users/jargrett/Desktop/Castilleja/Data Analysis/RMBL/2024 Argrett/Diversity 2023.csv", row.names=FALSE)


el.div.23 <- read.csv("Diversity 2023.csv")
el.div.23$year <- 2023
el.div.24 <- read.csv("Diversity 2024.csv")
el.div.24$year <- 2024

el.div.comb<- rbind.fill(el.div.23,el.div.24)
el.div.comb$year <- as.factor(el.div.comb$year)

#Run some models
div.lmm <- lmer(div ~ litter*removal + year + (1|block) + (1|pair), data = el.div.comb)
summary(div.lmm)
Anova(div.lmm)
emmip(div.lmm, removal ~ year)
emmeans(div.lmm, pairwise ~ removal|year)

rich.lmm <- lmer(rich ~ litter*removal + year + (1|block) + (1|pair), data = el.div.comb)
summary(rich.lmm)
Anova(rich.lmm)
emmip(rich.lmm, removal ~ year)
emmeans(rich.lmm, pairwise ~ removal|year)

even.lmm <- lmer(even ~ litter*removal + year + (1|block) + (1|pair), data = el.div.comb)
summary(even.lmm)
Anova(even.lmm)
emmip(even.lmm, removal ~ year)
emmeans(even.lmm, pairwise ~ year|removal)

#Standard error calculations
el.div <- el.div.comb %>% 
  group_by(removal, year) %>% 
  dplyr::summarise(mean= mean(div),
                   se = sd(div)/sqrt(n()))
el.rich <- el.div.comb %>% 
  group_by(removal, year) %>% 
  dplyr::summarise(mean= mean(rich),
                   se = sd(rich)/sqrt(n()))
el.even <- el.div.comb %>% 
  group_by(removal, year) %>% 
  dplyr::summarise(mean= mean(even),
                   se = sd(even)/sqrt(n()))

div.plot <- ggplot(data = el.div, aes(x = removal, y = mean, color = removal)) +
  stat_summary(fun=mean, colour="grey90", geom="line", aes(group = 1)) +
  geom_point(shape=18, size = 4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  facet_wrap(~year) + 
  scale_color_manual( values=c("#35978f", "#8c510a")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray", linewidth = 0.12)) +
  labs(x = "Castilleja", y = "Shannon diversity of co-occuring species") +
  theme(legend.position="none") +
  ylim(1.5,2.5)

div.plot

rich.plot <- ggplot(data = el.rich, aes(x = removal, y = mean, color = removal)) +
  stat_summary(fun=mean, colour="grey90", geom="line", aes(group = 1)) +
  geom_point(shape=18, size = 4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  facet_wrap(~year) + 
  scale_color_manual( values=c("#35978f", "#8c510a")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray", linewidth = 0.12)) +
  labs(x = "Castilleja", y = "Species Richness of co-occuring species") +
  theme(legend.position="none") +
  ylim(10,20)

rich.plot

even.plot <- ggplot(data = el.even, aes(x = removal, y = mean, color = removal)) +
  stat_summary(fun=mean, colour="grey90", geom="line", aes(group = 1)) +
  geom_point(shape=18, size = 4) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  facet_wrap(~year) + 
  scale_color_manual( values=c("#35978f", "#8c510a")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray", linewidth = 0.12)) +
  labs(x = "Castilleja", y = "Species Evenness of co-occuring species") +
  theme(legend.position="none") +
  ylim(0.5,1)

even.plot

diversity.plots <- ggarrange(div.plot, rich.plot, even.plot,
                             labels = c("A", "B","C"), 
                             nrow = 1)
diversity.plots  

alt.rich.plot <- ggplot(data = el.rich, aes(x = year, y = mean, color = removal)) +
  stat_summary(fun=mean, colour="grey90", geom="line", aes(group = 1)) +
  geom_point(shape=18, size = 4,
             position =  position_dodge(width = 0.5), width = 0.07) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),
                position =  position_dodge(width = 0.5), width = 0.07) +
  theme_pubr() +
  scale_color_manual( values=c("#35978f", "#8c510a")) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", 
                                    color = "gray", linewidth = 0.12)) +
  labs(x = "Sampling year", y = "Species Richness of co-occuring species") +
  theme(legend.position="none") +
  ylim(10,20)

alt.rich.plot

#--------------------------Multivariate analysis-------------------------------#
#we will now run a (Multivariate analysis)
#This allows us to look at the compositional differences between our sites,castilleja,etc.
#we are working towards Matrix format so we can take our castilleja matrix as our starting point
library(ggrepel)
library(vegan)
library(ggordiplots)



#First calculate distance matrix
dist <-vegdist(emerald.24.matrix, method="bray")


set.seed(20)
#Run NMDS on distance matrix
nmds <- metaMDS(dist, distance="bray", #use bray-curtis distance
                k=3, #2 dimensions
                try=500) #for publication I recommend 500)


nmds#stress value 0.14 which is below .2 so we need to investigate

ordiplot(nmds, type="text", display="sites")

nmds.scores <- as.data.frame(vegan::scores(nmds))

NMDS <- cbind(el.div.24,nmds.scores) #final dataset


perm <- adonis2(dist ~ litter*removal, data = NMDS, permutations=9999)
perm


ggplot(NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=litter , shape=removal), size = 3, alpha = 0.8) +
  scale_color_manual( values=c("#B308F6", "#35978f","#E1BE6A","#000000")) +
  coord_equal() +
  theme_bw()
