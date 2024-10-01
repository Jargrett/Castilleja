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

summary(gbif_cl)
occurance_ga <- subset(gbif_cl, stateProvince == 'Georgia')
write.csv(occurance_ga, "/Users/jargrett/Desktop/Castilleja/Data Analysis/Bot RA/Georgia Occurances.csv", row.names=FALSE)
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



