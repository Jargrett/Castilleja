#Biocrust example code
#AC 2100127

#Set up---------------------------
#Set working directory, clean environment, etc.
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Anny Stats/Lab 4")

#Load packages
library(car)
library(ggpubr)
library(emmeans)
library(rstatix)

#load data
biocrust <- read.csv("biocrust.csv")

#Understand the data----------
str(biocrust)
#Plot: 40 plots total, 1-20 are in site C, 21-40 are in site G
#Location: g=grass, s=shrub, i=interspace; collected within plots
#Chlor: chlorophyll concentration in surface soil
#Scyt: scytonemin concentration in surface soil
#Treatment: control vs. stomp applied to each plot
#Site: C=creosote site, G=grassland site

table(biocrust$Treatment,biocrust$Location,biocrust$Site)
#Unbalanced design even without considering NAs in the response
#There are no shrubs (s) in the grassland site (G)
#Twice as many interspace (i) samples as samples associated with grasses (g) or shrubs (s)

#Use site=C and response=chlor as example----------------------
#Subset data
data.C<-subset(biocrust,Site=="C")

#Prelim viz
hist(data.C$Chlor_ugpergsoil)#Very not normal, but let's reserve judgement
boxplot(data.C$Chlor_ugpergsoil~data.C$Location*data.C$Treatment)#Some outliers in g.control

#Does chlorophyll concentration differ among locations and treatments?

#Try raw data
m1<-lm(Chlor_ugpergsoil~Location*Treatment, data = data.C)
hist(m1$residuals)#not very normal
qqPlot(m1)#outliers

#log transform
#check for zeros
which(data.C$Chlor_ugpergsoil==0)#zeros exist
#calculate minimum non-zero value to add to all data except NAs to avoid log(0)
min.nz <- min(data.C$Chlor_ugpergsoil[which(data.C$Chlor_ugpergsoil>0)])
#model with log-transformed data
m2<-lm(log(Chlor_ugpergsoil+min.nz)~Location*Treatment, data=data.C)
hist(m2$residuals)
qqPlot(m2)#much better
plot(m2)#seems okay

#do another model with different effects order for comparison
m3<-lm(log(Chlor_ugpergsoil+min.nz)~Treatment*Location, data=data.C)

#Use m2 and m3 to explore how SS affects results-------------
#Summary of lm not exactly helpful
summary(m2) #F-statistic: 10.71 on 5 and 152 DF,  p-value: 7.791e-09 (note 2 NAs)
summary(m3) #F-statistic: 10.71 on 5 and 152 DF,  p-value: 7.791e-09 (same)

#What is the "intercept" here?
head(model.matrix(m2))
head(data.C[,c("Location","Treatment")])
#g-control

#Compare Type I: SS and results change depending on order
summary(aov(log(Chlor_ugpergsoil+min.nz)~Location*Treatment, data=data.C)) #Same as anova(m2)
#Location: F=20.647 P=1.17e-08
#Treatment: F=9.220  P=0.00282
#Interaction: F=1.526  P=0.22075
summary(aov(log(Chlor_ugpergsoil+min.nz)~Treatment*Location, data=data.C)) #Same as anova(m3)
#Location: F=20.374 P=1.45e-08
#Treatment: F=9.766  P=0.00213 
#Interaction: F=1.526  P=0.22075

#Type 2 does not depend on effects order
Anova(m2, type=2)
#Location: F=20.3742 P=1.448e-08
#Treatment: F=9.2199  P=0.00282
#Interaction: F=1.5258  P=0.220749
Anova(m3, type=2) 
#Location: F=20.3742 P=1.448e-08
#Treatment: F=9.2199  P=0.002818
#Interaction: F=1.5258  P=0.220749

#Go back to ppt to talk about reporting

#Visualization---------------------
#Plot means and SEs at original untransformed scale
p<-ggerrorplot(data=data.C, x="Location", y="Chlor_ugpergsoil",
            color="Treatment",palette=c("darkgray","black"),
            desc_stat = "mean_se",
            position = position_dodge(0.3),
            xlab="substrate type", ylab="chlorophyll concentration (ug/g)",
            size=0.7
            )

p
#Explore issues with plotting transformed data vs. untransformed data
# add = "jitter", add.params = list(color = "darkgray"),
#yscale="log2",

#Post hoc comparison of means---------------------
#We can technically still do Tukey HSD, but remember it only works on an aov object
#And aov() conducts Type I ANOVA, which is inappropriate for this dataset
#Try anyway
TukeyHSD(aov(log(Chlor_ugpergsoil+min.nz)~Location*Treatment, data=data.C))
#Lots of output, does it all relate to a valid hypothesis?

#A more flexible option: emmeans package (estimated marginal means)
#Example similar to Tukey above
emmeans(m2, specs=pairwise~Treatment*Location, adjust="tukey")

#Q: Does substrate type change stomp effect?
em<-emmeans(m2, specs=pairwise~Treatment|Location)
em #it does not appear that p values are adjusted (because emmeans defaults "family to each group)
em$contrasts %>%
  rbind(adjust="tukey")#default is bonferroni, very conservative
#Yes! Stomping only decreases chlorophyll in the grass and interspace substrates, but not near shrubs

#It might be nice to have the estimated means on the response (not logged) scale
em2<-emmeans(m2, specs=pairwise~Treatment|Location, type="response")
em2$contrasts %>%
  rbind(adjust="tukey")
#test statistics are the same, the estimated means are now on the original scale

#We can add these P values to our plot using the rstatix and ggpubr packages
#Add to original plot using geom_bracket()
p+geom_bracket(
  xmin = c(1,2,3)-0.07, xmax = c(1,2,3)+0.07, y.position = c(0.5, 0.2, 0.3),
  label = c("p = 0.03", "p = 0.04", "p = 0.99"),
  tip.length = c(0.01, 0.01), vjust = -0.4
) 

#Global model to demonstrate it doesn't work--------------------------------------
m0<-lm(Chlor_ugpergsoil~Location*Treatment*Site, data = biocrust)
qqPlot(m0)
min(biocrust$Chlor_ugpergsoil[which(biocrust$Chlor_ugpergsoil>0)])
m01<-lm(log(Chlor_ugpergsoil+0.01346004)~Location*Treatment*Site, data = biocrust)
qqPlot(m01)
Anova(m01)
#Error "model has aliased coefficients"
#This means that there are levels of variables that only occur for levels of other variables (e.g. not full factorial)
#In our case, location=s only occurs at site = C, so the effect of shrub and shrubland site are confounded
#Therefore, must analyze sites separately, or only consider grass vs. interspace and ignore shrubs

