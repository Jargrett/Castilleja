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
# here we will be working with the count column
# we will need to conververt data to a matrix format
emerald.24 <- EL.comb%>% filter (year == "2024")
emerald.24 <- emerald.24[!(emerald.24$functional_group %in% "environmental"),]
emerald.24 <- emerald.24[!(emerald.24$code %in% "CASE"),]
emerald.24.cov <- subset(emerald.24, select = c('plot','code','cover'))

emerald.23 <- EL.comb%>% filter (year == "2023")
emerald.23 <- emerald.24[!(emerald.23$functional_group %in% "environmental"),]
emerald.23 <- emerald.24[!(emerald.23$code %in% "CASE"),]
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
div.lmm <- lmer(div ~ removal*year + (1|block) + (1|pair), data = el.div.comb)
summary(div.lmm)
Anova(div.lmm)
emmip(div.lmm, removal ~ year)
emmeans(div.lmm, pairwise ~ removal|year)

rich.lmm <- lmer(rich ~ removal*year + (1|block) + (1|pair), data = el.div.comb)
summary(rich.lmm)
Anova(rich.lmm)
emmip(rich.lmm, removal ~ year)
emmeans(rich.lmm, pairwise ~ removal|year)

even.lmm <- lmer(even ~ removal*year + (1|block) + (1|pair), data = el.div.comb)
summary(even.lmm)
Anova(even.lmm)
emmip(even.lmm, removal ~ year)
emmeans(even.lmm, pairwise ~ year|removal)

diversity.plot <- ggplot(el.div.comb, aes(x = removal, y = div)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (removal), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~year) +
  scale_color_manual(values=c( "khaki3", "burlywood4")) +
  labs(x = "Castilleja septentrionalis", y = "Shannon Diversity") +
  ylim(1,3)

diversity.plot

richness.plot <- ggplot(el.div.comb, aes(x = removal, y = rich)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (removal), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~year) +
  scale_color_manual(values=c( "khaki3", "burlywood4")) +
  labs(x = "Castilleja septentrionalis", y = "Species Richness") +
  ylim(5,22)

richness.plot

evenness.plot <- ggplot(el.div.comb, aes(x = removal, y = even)) +
  stat_summary(aes(group = pair), geom = "line", fun.y = mean, col ="ivory3") +
  geom_point(aes(color = (removal), size = 1, alpha = 2), show.legend = FALSE) + 
  stat_summary(fun=mean, geom = "crossbar", position = position_dodge(1), linewidth = 1, width = 0.25, col = "grey34") +
  theme_pubr() +
  facet_wrap(~year) +
  scale_color_manual(values=c( "khaki3", "burlywood4")) +
  labs(x = "Castilleja septentrionalis", y = "Species Evenness") +
  ylim(0.5,1)

evenness.plot


div.plot <- ggarrange(diversity.plot, richness.plot, evenness.plot,
                                  labels = c("A", "B","C"), 
                                  nrow = 1, common.legend = TRUE, legend = "bottom")
div.plot

#--------Composition Analysis (Multivariate analysis) -------#
library(ggrepel)
# non-metric multidimensional scaling (NMDS)
# for this we need to take a matrix so we will retake code from before
emerald.24 <- EL.comb%>% filter (year == "2024")
emerald.24 <- emerald.24[!(emerald.24$functional_group %in% "environmental"),] #remove rock,bare, and Litter

emerald.matrix <- matrify(emerald.24.cov)
set.seed(20)#This sets the random start seed so that we are guaranteed to all get the same outputs

#First calculate distance matrix
emerald.dist <-vegdist(emerald.matrix, method="bray")

#Run NMDS on distance matrix
emerald.nmds<-metaMDS(emerald.dist,distance="bray", #use bray-curtis distance
                   k=3, #2 dimensions
                   try=500 #force it to try 100 different times (default is 20, for publication I recommend 500)
)
emerald.nmds #Gives stress value of best solution, aim for <0.2
#Check the fit

stressplot(emerald.nmds) #Also called a "Shepard diagram"
#We want this plot to show a monotonic relationship

ordiplot(emerald.nmds, type="text", display="sites")

#Can output ordination coordinates and use ggplot for nicer graphics
nmds.scores <- as.data.frame(vegan::scores(emerald.nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
#Add to mite.info dataframe
#combined data set with Plot Data File and calculated values
el.data <- read.csv("Diversity 2024.csv") #importing metadata
el.NMDS <- cbind(el.data,nmds.scores) #final dataset

#Build NMDS figure in ggplot
ggplot(el.NMDS, aes(NMDS1, NMDS2)) +
  geom_point() +
  geom_text_repel(aes(NMDS1, NMDS2, label = pair), size = 3)#avoids text overlapping

NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 6), replicate(6, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 6),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS Scree plot")
  for (i in 1:6) {
    points(rep(i + 1,6),replicate(6, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}
#looks like k=3 is the best for our data
NMDS.scree(emerald.dist)

#PERMANOVA to test whether composition differs by shrub cover or topography---------
adonis2(emerald.dist~removal*block*litter, data = el.NMDS, permutations=999)


ggplot(el.NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(color= litter, shape=removal)) +
  coord_equal() +
  theme_bw()

