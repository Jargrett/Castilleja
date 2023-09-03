#Field Plant Soil Feedback Experiment mixed model example
#AC 210310


#Set up---------------------------
#Set working directory, clean environment, etc.
setwd("C:/Users/jargr/Dropbox/PC/Desktop/Anny Stats/PSF Example 3_22")
#load libraries
library(car)
library(lme4)
library(lmerTest)
library(performance)
library(see)
library(ggpubr)
library(emmeans)

#load data
BOGR.dat <- read.csv("SEV_PSF_BOGR_data.csv")

#Prelim visualizations--------------------
hist(BOGR.dat$logMass.harv)
boxplot(BOGR.dat$logMass.harv~BOGR.dat$PSF*BOGR.dat$State)

#Visualize the effect of repeated measures
#Use ggpubr package for quick plotting
ggscatter(data=BOGR.dat, x="Census", y="logMass.harv", color = "PSF", shape="State",
          facet.by="ID")

#Visualize spatial blocking effects
ggscatter(data=BOGR.dat, x="Census", y="logMass.harv", color = "PSF", shape="State",
          facet.by="Block")

#Simple mixed model for repeated sampling and spatial blocking-------------
m1.BOGR<-lmer(logMass.harv~PSF*State+Census+(1|Block)+(1|ID),data=BOGR.dat)
check_model(m1.BOGR)
summary(m1.BOGR)
Anova(m1.BOGR)

#Quick results plot
#Exponentiate (back-transform) logMass.harv to plot biomass at observed scale
BOGR.dat$Mass.harv<-exp(BOGR.dat$logMass.harv)
ggerrorplot(data=BOGR.dat, x="State", y="Mass.harv", color="PSF",
            ylab = "Aboveground biomass (g)")
#These SE's are directly calculated from data, thus extremely small since it assumes all 445 observations are independent

#Post-hoc comparisons
emmeans(m1.BOGR, pairwise~PSF|State, adjust="Tukey")#Specifically only compare home vs. away within each State

#Conclude that the effect of plant-soil feedback environment on BOGR growth depends on the historical dynamics

#Alternative models--------------------------------------
#Compare to the nested version if I had nested ID's instead
m1A.BOGR<-lmer(logMass.harv~PSF*State+Census+(1|Block/nested.ID),data=BOGR.dat)
summary(m1A.BOGR)#Gives same results as above


#We can also try this using Census as a random effect instead, residuals a little worse, results similar
m2.BOGR<-lmer(logMass.harv~PSF*State+(1|Census)+(1|Block)+(1|ID),data=BOGR.dat)
check_model(m2.BOGR)
summary(m2.BOGR)
Anova(m2.BOGR)

#We can also try introducing random slopes, residuals a little worse, results similar
m3.BOGR<-lmer(logMass.harv~PSF*State+Census+(State|Block)+(Census|ID),data=BOGR.dat)
check_model(m3.BOGR)
summary(m3.BOGR)
Anova(m3.BOGR)

#Quick generalized linear mixed model GLMM example------------
#Analyze tiller number (technically count data)
m4.BOGR<-glmer(Tiller~PSF*State+Census++(1|Block)+(1|ID),data=BOGR.dat, family = poisson)
check_model(m4.BOGR)
#Check for overdispersion 
resid.ssq <- sum(residuals(m4.BOGR,type="pearson")^2)  
resid.df <- nrow(subset(BOGR.dat,!is.na(Tiller) ))-length(coef(m4.BOGR)) 
resid.ssq/resid.df #not bad
#Check results
summary(m4.BOGR)
Anova(m4.BOGR)

#Quick autocorrelation example---------------------
library(nlme)
m5.BOGR<-lme(logMass.harv~PSF*State+Census, random = list(~1|Block,~1|ID),data=BOGR.dat)
summary(m5.BOGR)#Note that lme() cannot do crossed random effects right
plot(ACF(m5.BOGR, resType = "normalized"), alpha = 0.05)
#ACF assumes that your time points (censuses here) are equidistant

#Add AR1 correlation structure
m5.AR.BOGR<-lme(logMass.harv~PSF*State+Census, random = list(~1|Block,~1|ID),correlation=corAR1(form=~Census),data=BOGR.dat)
summary(m5.AR.BOGR)#Lower AIC value compared to before
