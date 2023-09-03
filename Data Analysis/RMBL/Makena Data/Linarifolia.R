# Deer Creek linariifolia analysis
# Initiated: 7/20/23
# Completed: TBD

#####################
#                   # - Impact of Castilleja presence on:
#      Question     # - richness and evenness (Shannon Diversity Index)
#                   # - 
##################### - Also would be good too look at: Individual species cover

# Set working directory (Workspace)
# This can be done manually or through session -> set working directory
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Data Analysis/RMBL/Makena Data")

# Load-in packages
# I will explain these if needed
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(tidyverse) # for data working
library(vegan) # for diversity
library(ggplot2) # plotting
library(ggpubr) # plotting
library(car)
library(emmeans)
library(rstatix)
library(ggrepel)
library(devtools)
library(pairwiseAdonis)
library(indicspecies)
library(simboot)

# importing dataframe

dc <- read.csv("Combined deer creek data - Individuals.csv")
dc.counts <- dc[ -c(1:5,27,28)] # removing unnecessary columns 
dc.env <- subset(dc, select=c(2,4,5)) # holding our environmental stuff
dc.covers <- read.csv("Combined deer creek data - Cover.csv")
dc.covers <- dc.covers[ -c(1:6,27,28)]
#--------------------Diversity Analysis--------------------#
# calculating diversity metrics for our linariifolia sites

# Calculating the species richness for plots
rich <- specnumber(dc.counts)

# Calculating Shannon diversity for plots
div <- diversity(dc.counts, index = "shannon")
div.cov <- diversity(dc.covers, index = "shannon")

# Calculating species evenness for plots
even <- diversity(dc.counts, index = "shannon") / log(specnumber(dc.counts))     
even.cov <- diversity(dc.covers, index = "shannon") / log(specnumber(dc.counts))     
#combined data set with environmental and calculated values
dc.div <- cbind(dc.env,rich,div,even)
dc.cov.div <- cbind(dc.env,div.cov,even.cov,rich)

#export this dataset for combined analysis
write.csv(dc.cov.div, "C:\\Users\\jargr\\Dropbox\\PC\\Desktop\\Data Analysis\\RMBL\\Makena Data\\CALI Diversity.csv")

even.lm <- lm(even.cov ~ Treatment*Site, data = dc.cov.div)
summary(even.lm)
plot(even.lm)
Anova(even.lm)
hist(even.lm$residuals)

mean(dc.cov.div$div.cov)


#fit a GLM
rich.glm <- glm(rich ~ Treatment*Site, family= "poisson", data = dc.div)

hist(rich.glm$residuals) # pretty normal
plot(rich.glm)# seems okay

# We can get a summary of the model:
summary(rich.glm) 
Anova(rich.glm)

#fit a GLM


div.lm <- lm(div.cov ~ Treatment*Site, data = dc.cov.div)
summary(div.lm)
plot(div.lm)
Anova(div.lm)
h <- aov(div.cov ~ Treatment, data = dc.cov.div)
summary(h)

rich.em<-emmeans(rich.glm, specs=pairwise~Treatment, type="response")
rich.em$contrasts %>%
  rbind(adjust="Sidak")

rich.em

#----------------Visuals!!!----------------#

# First richness:
rich.plot <- ggplot(dc.div, aes(x = Treatment, y = rich, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Treatment", y = "Species Richness") +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3) +
  theme_classic2() +
  geom_text(x=1.5,y=14.5,label="**") +
  ylim(2.5,15) +
  geom_segment(aes(x=1,xend=2,y=14,yend=14))
  
rich.plot

div.plot <- ggplot(dc.cov.div, aes(x = Treatment, y = div.cov, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Treatment", y = "Shannon") +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3) +
  theme_classic2()+
  ylim(0.5,3)
div.plot

site.even.plot <- ggplot(dc.cov.div, aes(x = Site, y = even, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Population", y = "Species eveness") +
  geom_jitter(shape=16, alpha = 0.3) +
  theme_classic2()

site.even.plot

#---------------------Cover Analysis---------------------#
# We will use a PERMANOVA (adonis2() in the vegan package)
# to do this we need to convert our matrix values to a distance matrix
# NMDS First
# NMDS is a 2 dimensional representation of our data that show similarities dissimilarities between our data
dc.cover <- read.csv("Combined deer creek data - Cover.csv")
dc.cover[is.na(dc.cover)] = 0
dc.cover2 <- dc.cover[ -c(1:3,5)]
dc.cover.matrix <- dc.cover2 %>% remove_rownames %>% column_to_rownames(var="Plot")

dc.cover.info <- subset(dc.cover, select=c(1,2,4,5))
dc.cover.info$Treatment <- as.factor(dc.cover.info$Treatment)
is.na(dc.cover.matrix)
set.seed(20)

dc.cover.dist<-vegdist(dc.cover.matrix, method="bray")

dc.cover.nmds<-metaMDS(dc.cover.dist,distance="bray", k=2,try=500)
dc.cover.nmds
stressplot(dc.cover.nmds)

#PERMANOVA
perm.cover <- adonis2(dc.cover.dist~Treatment*Site, data = dc.cover.info, permutations=999)
perm.cover

pairwise.adonis2(dc.cover.dist~Treatment*Site, data = dc.cover.info, permutations=999)

dc.nmds.scores <- as.data.frame(vegan::scores(dc.cover.nmds))
dc.cover.info<-cbind.data.frame(dc.cover.info,dc.nmds.scores) 

ggplot(dc.cover.info, aes(NMDS1, NMDS2, color= Treatment, shape= Site, linetype = Site)) +
  geom_point() +
  stat_ellipse() 

# PERMANOVA RESULTS VISUALIZATIONS
# we first need to look at some of the raw data
# lets convert our data to list format
dc_long <- dc.cover %>% pivot_longer(cols='Bare.ground':'Wyethia.amplexicaulis',names_to="Taxa",values_to="Cover")
dc_long2 <- dc %>% pivot_longer(cols='Agastache.urticifolia':'Wyethia.amplexicaulis',names_to="Taxa",values_to="Count")

# subset data to values where cover !=0
# we will use this to run linear models on target functional group by treatment
dc.cover.long <- filter(dc_long, Cover > 0)
dc.count.long <- filter(dc_long2, Count > 0)

# merge this data with functional group and count data 
dc.species <- read.csv("Deer creek species list.csv")
dc.functional <- dc.species[ -c(1:3,7)]
dc.cover.taxa <- merge(dc.cover.long,dc.functional,by="Taxa")
dc.count.taxa <- merge(dc.count.long,dc.functional,by="Taxa")
dc.count.taxa2 <- subset(dc.count.taxa,select=c(1,7))
dc.taxa <- merge(dc.cover.taxa,dc.count.taxa2,by="Taxa")

# now to plot
# lets plot functional group by treatment
dc.taxa$Taxa <- as.factor(dc.taxa$Taxa)
dc.taxa$Functional.group <- as.factor(dc.taxa$Functional.group)
dc.taxa$Life.history <- as.factor(dc.taxa$Life.history)


treatment.by.fgroup <- ggplot(dc.taxa, aes(x = Treatment, fill = Functional.group, y = Cover)) + 
  geom_bar(stat = "identity") +
  labs(x = "Treatment", y = "Total Percent Cover") +
  facet_wrap(~Functional.group)
treatment.by.fgroup

#---------------------Species Cover Analysis---------------------#
# First we will do some preliminary linear regressions

m <- lm(Cover~Treatment*Site, data = subset(dc.cover.taxa, Functional.group ==""))
summary(m)
Anova(m)

b <- lm(Cover~Treatment*Site, data = dc.cover.taxa)
Anova(b)
summary(b)

plot.species.cover<-ggplot(subset(dc.taxa,Functional.group %in% "Grass"), aes(x=Treatment, y= Cover)) +
  geom_boxplot()
plot.species.cover

plot.species.cover<-ggplot(dc.cover.taxa, aes(x=Site, y= Cover, fill = Treatment )) +
  geom_boxplot() +
  ylim(0,.2)
plot.species.cover
