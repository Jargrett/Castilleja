#Species codes
#initiated: 9/26/24
setwd("/Users/jargrett/Desktop/Castilleja/Data Analysis/Bot RA")
#install.packages("rgbif")
library(rgbif)
#install.packages("countrycode")
library(countrycode)
#install.packages("CoordinateCleaner")
library(CoordinateCleaner)
#install.packages("sf")
library(sf)
library(tidyverse)
user <- "jargrett" # your gbif.org username 
pwd <- "Korra13506!" # your gbif.org password
email <- "jordan.argrett@uga.edu" # your email 

gbif_data <- read.csv("GBIF1.csv")

gbif_clean <- gbif_data %>% select(species,decimalLongitude, 
                decimalLatitude, countryCode, stateProvince, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters,
                year, basisOfRecord, institutionCode, datasetName)
summary(gbif_clean)

gbif_cleaner <- gbif_clean %>%
  filter(!is.na(decimalLongitude)) %>%
  filter(!is.na(decimalLatitude))

wm <- borders("world", colour = "gray50", fill = "gray50")
ggplot() +
  coord_fixed() +
  wm +
  geom_point(data = gbif_cleaner,
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred",
             size = 0.5) +
  theme_bw()

gbif_cleaner$countryCode <-  countrycode(gbif_cleaner$countryCode, 
                                origin =  'iso2c',
                                destination = 'iso3c')

gbif_cleaner <- data.frame(gbif_cleaner)
flags <- clean_coordinates(x = gbif_cleaner, 
                           lon = "decimalLongitude", 
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("capitals", "centroids",
                                     "equal", "zeros", "countries"))
summary(flags)
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")

gbif_cl <- gbif_cleaner[flags$.summary,]
write.csv(gbif_cl, "/Users/jargrett/Desktop/Castilleja/Data Analysis/Bot RA/Total Occurances.csv", row.names=FALSE)

summary(gbif_cl)
occurance_us <- subset(gbif_cl, countryCode == 'USA')
write.csv(occurance_us, "/Users/jargrett/Desktop/Castilleja/Data Analysis/Bot RA/US Occurances.csv", row.names=FALSE)
occurance_ga <- subset(gbif_cl, stateProvince == 'Georgia')
write.csv(occurance_ga, "/Users/jargrett/Desktop/Castilleja/Data Analysis/Bot RA/Georgia Occurances.csv", row.names=FALSE)
cm <- borders("state", colour = "gray50", fill = "gray50")
ggplot() +
  coord_fixed() +
  cm +
  geom_point(data = occurance_ga,
             aes(x = decimalLongitude, y = decimalLatitude, fill = basisOfRecord),
             colour = "darkred",
             size = 0.5) +
  theme_bw()
#Species codes:
#GBIF.org (01 October 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.2f8w82
#Parthenium integrifolium
#Silphium compositum
#Packera anonyma
#Asclepias tuberosa
#Coreopsis major
#Schizachyrium scoparium
#Eryngium yuccifolium

#taxon keys
PAIN <- name_backbone("Parthenium integrifolium")$usageKey
SICO <- name_backbone("Silphium compositum")$usageKey
PAAN <- name_backbone("Packera anonyma")$usageKey
ASTU <- name_backbone("Asclepias tuberosa")$usageKey
COMA <- name_backbone("Coreopsis major")$usageKey
SCSC <- name_backbone("Schizachyrium scoparium")$usageKey
ERYU <- name_backbone("Eryngium yuccifolium")$usageKey

species <- data.frame(
  taxonKey= c(PAIN,SICO,PAAN,ASTU,COMA,SCSC,ERYU))

occ_download(pred("taxonKey", 3086802),format = "SIMPLE_CSV", user=user,pwd=pwd,email=email)
occ_download_wait('0031966-240906103802322')
PAIN <- occ_download_get('0031967-240906103802322') %>%
  occ_download_import()

#---------#
### install all the packages we need
# Note: you need R version 4.0 or higher for dismo
library(raster)
#install.packages("megaSDM")
library(sp)
library(megaSDM)
library(sdm)
library(geodata)
library(ENMTools)
library(reshape2)
#devtools::install_github("brshipley/megaSDM", build_vignettes = TRUE)
### load environmental data from worldclim

# first create a folder in your working directory called "data"
# you'll download the environmental data into this folder

# note: you need to set download=TRUE the first time you download these files
# but you can change to download=F for any subsequent time you re-fun the script

# we are downloading bioclim variables at a 2.5 degree minute resolution
bio_curr <- worldclim_global(var='bio', res=0.5, download=TRUE, path="data/")

# we are downloading forecasted values of the bioclim variables
# from SSP585 for the time frame 2061-2080 and the model ACCESS-CM2
# (in practice you would look at a range of SSP scenarios and models)
bio_fut <- cmip6_world(var='bioc', res=2.5, 
                       ssp=585, time='2061-2080', model='ACCESS-CM2', download=TRUE, path="data/")

# look at the options you have for data here:
?worldclim_global()
?cmip6_world()


### read in data file of downsampled (200 km) GBIF occurrences
### for Chamaecrista fasciculata

occ<-read.csv("Georgia Occurances.csv")
head(occ) # coordinates of observed occurrences in N. America
summary(occ)

# establish spatial extent of the study
max.lat <- ceiling(max(occ$decimalLatitude))
min.lat <- floor(min(occ$decimalLatitude))
max.lon <- ceiling(max(occ$decimalLongitude))
min.lon <- floor(min(occ$decimalLongitude))

# create a map of the occurrences
world <- map_data("world")

ggplot()+geom_polygon(data=world,aes(x=long,y=lat,group=group),fill='gray90',color='black')+
  coord_fixed(xlim=c(min.lon, max.lon),ylim=c(min.lat, max.lat))+
  geom_point(data=occ,aes(x=decimalLongitude,y=decimalLatitude),color='olivedrab')+
  theme_bw()+ylab('Latitude')+xlab('Longitude')

### generate background points (pseudo-absence data)

# crop climate data to the spatial extent of the study
study <- extent(min.lon, max.lon, min.lat, max.lat)
bio_curr2 <- crop(bio_curr, study)
bio_fut2 <- crop(bio_fut, study)

# here we are sampling points in different cells than our presences 
# and using the Varela method for environmental subsampling
BackgroundPoints(spplist=c('Cfasc'),envdata=bio_curr2, output='data', nbg=nrow(occ), method='Varela')
bg.points <- read.csv('data/Cfasc_background.csv')
head(bg.points)

# we only need to keep the lat/longs
bg.points <- as.data.frame(bg.points)[,c('x','y')]
colnames(bg.points) <- c('decimalLongitude', 'decimalLatitude')
head(bg.points)

# visualize background points
ggplot()+geom_polygon(data=world,aes(x=long,y=lat,group=group),fill='gray90',color='black')+
  coord_fixed(xlim=c(min.lon, max.lon),ylim=c(min.lat, max.lat))+
  geom_point(data=occ,aes(x=decimalLongitude,y=decimalLatitude),color='olivedrab')+
  geom_point(data=bg.points,aes(x=decimalLongitude,y=decimalLatitude),color='gray30',pch=1)+
  theme_bw()+ylab('Latitude')+xlab('Longitude')

# now combine into one dataframe with presence/absence as a binary variable
occ$presence <- 1
bg.points$presence <- 0
all_points <- rbind(occ, bg.points)
head(all_points)

#### Now get data into the formats we need for sdms

# re-crop climate layters to the extent of all_points
max.lat2 <- ceiling(max(all_points$latitude))
min.lat2 <- floor(min(all_points$latitude))
max.lon2 <- ceiling(max(all_points$longitude))
min.lon2 <- floor(min(all_points$longitude))

study2 <- extent(min.lon2, max.lon2, min.lat2, max.lat2)
bio_curr3 <- crop(bio_curr, study2)
bio_fut3 <- crop(bio_fut, study2)

# give the climate variables the same names
names(bio_curr3)
names(bio_fut3)
names(bio_fut3)=names(bio_curr3)=paste('bio',1:19,sep='') # note: this works because variables 1-19 are in the same order in both

# Now we need to turn our climate rasters into stacks
bio_curr3 <- stack(bio_curr3)
bio_fut3 <- stack(bio_fut3)

# convert to SpatialPointsDataFrame
species <- SpatialPointsDataFrame(coords = all_points[, c("longitude", "latitude")],
                                  data = all_points,
                                  proj4string = CRS("+proj=longlat +datum=WGS84"))

e <- as(extent(bio_curr3), 'SpatialPolygons')
crs(e) <- crs(species)

# visualize current and future environmental data
plot(bio_curr3)
plot(bio_fut3)

# look at correlations among climate variables
raster.cor.plot(bio_curr3)
raster.cor.plot(bio_fut3)

# these are pretty correlated -- we might want to think about using a subset of predictors

###### create the data object that you'll pass to the sdm function
dat <- sdmData(train = species, predictors = bio_curr3)

############# Fit some distribution models!

### first we have to decide which climate variables to use as predictors
# here we're using a subset we think would be important biologically
formula_sub <- presence ~ bio3 + bio9 + bio10 + bio17 + bio18

# here we're using all 19
formula_all <- presence ~ bio1+bio2+bio3+bio4+bio5+bio6+bio7+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19

# you can make whatever formula you want by just listing the bio variables you want to include

### next we have to decide how to split data into training vs. testing
# there are lots of ways to evaluate fit. Here we're subsampling 30% of the data
# to use as testing data, and repeating that 5 times
# in practice you would likely use more reps, but it takes a while to run
percent <- 30
rep <- 5

### next we have to decide what algorithm to use
# options are gam, glm, gbm, svm, rf, tree, mars, mda, fda, brt, maxent

# let's start with random forest using all 19 predictors
mod.rf <- sdm(formula = formula_all, data = dat, methods = 'rf', 
              replication = "sub", test.percent = percent, n=rep)


# look at the model info
mod.rf

# look at the predictions based on current climate
# we are taking the mean across all the different runs
mod.rf.predict <- predict(mod.rf, bio_curr3, ext=study2, mean=T)
plot(mod.rf.predict, main='Current climate')

# look at the relative importance of each climate variable
# which are important for the biology of the species?
vi <- getVarImp(mod.rf, id=1:rep)
plot(vi) 

# look at the partial response curves to each climate variable
# (holding others constant at their mean)
# do these make sense biologically?
rcurve(mod.rf, id=1:rep, mean = T, confidence = T)

# look at the fit of the model
getEvaluation(mod.rf) 

# look at the predictions based on future climate
mod.rf.forecast <- predict(mod.rf, newdata=bio_fut3, ext=study, mean=T)
plot(mod.rf.forecast, main='SSP 585')

############ That was with one algorithm, but in practice we would want to do many

# make an ensemble model with many algorithms
# note: you can change the methods argument to include whichever algorithms you want
mod.all<-sdm(formula = formula_all, data = dat, methods = c("gam", "brt", "mars", "svm", "glm", "rf"), replication = 'sub', test.percent = percent, n=rep)
mod.all

#### Now we need to combine our predictions in some way
# there are many ways to do this
# here, we are taking a weighted average across models based on each model's TSS

# create an ensemble prediction based on current climate
mod.all.ensemble <- ensemble(mod.all,
                             newdata = bio_curr3,
                             setting = list(method ='weighted', stat = 'TSS', opt = 2),
                             overwrite = TRUE)

# create an ensemble prediction based on future climate
mod.all.ensemble.future <- ensemble(mod.all,
                                    newdata = bio_fut3,
                                    setting = list(method ='weighted', stat = 'TSS', opt = 2),
                                    overwrite = TRUE)

plot(mod.all.ensemble, main = "Bioclim Contemporary Ensemble")

plot(mod.all.ensemble.future, main = "SSP245 2061-2080 Ensemble")

#### Sometimes we want to convert continuous suitability scores into binary presence/absence maps
# here we are taking the average presence/absence of each model (fraction of models that predict presence)
mod.all.ensemble.thresh <- ensemble(mod.all,
                                    newdata = bio_curr3,
                                    setting = list(method ='pa', stat = 'TSS', opt = 2),
                                    overwrite = TRUE)

# create an ensemble prediction based on future climate
mod.all.ensemble.future.thresh <- ensemble(mod.all,
                                           newdata = bio_fut3,
                                           setting = list(method ='pa', stat = 'TSS', opt = 2),
                                           overwrite = TRUE)

plot(mod.all.ensemble.thresh, main = "Bioclim Contemporary Ensemble Presence")

plot(mod.all.ensemble.future.thresh, main = "SSP585 2061-2080 Ensemble Presence")


