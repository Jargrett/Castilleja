# Download_Occurrence_Data_Updated.R
# ---------------------------------------------------------------
# Purpose: Download and preview species occurrence records using iDigBio and gatoRs.
# Created by ML Gaynor
# ---------------------------------------------------------------
setwd("~/Desktop/Castilleja/Data Analysis/Accession work")
## ---- Load Required Packages ----
library(ridigbio)     # Interface to iDigBio API
library(gatoRs)       # Unified taxonomic/occurrence tools
library(leaflet)      # Interactive mapping

## ---- A) Download from iDigBio ----

# Search for specific species (Galax urceolata) 
iDigBio_GU <- idig_search_records(rq = list(scientificname = "Macranthera flammea")) #rq = research query

# Search for all Diapensiaceae occurrences (limit to 1000 for example)
iDigBio_GU_family <- idig_search_records(rq = list(family = "Orobanchaceae"), limit = 1000)


# Search with a geographic bounding box (e.g., eastern USA extent)
rq_input <- list(
  scientificname = list(type = "exists"),
  family = "Orobanchaceae",
  geopoint = list(
    type = "geo_bounding_box",
    top_left = list(lat = 37.99616267972814, lon = -97.73437500000001),
    bottom_right = list(lat = 25.799891182088334, lon = -71.36718750000001)
  )
)
iDigBio_GU_family_USA <- idig_search_records(rq_input, limit = 1000)

# Save iDigBio results as CSV
write.csv(iDigBio_GU, "iDigBio_MF_2025_06_27.csv", row.names = FALSE) #Name csv pull with source_species_date
write.csv(iDigBio_GU_family, "iDigBio_GU_Oro_2025_06_27.csv", row.names = FALSE)

## ---- B) Download Using gatoRs ----

# Define synonym lists for each focal species
Shortia_galacifolia <- c("Macranthera flammea", "Sherwoodia galacifolia")
Galax_urceolata <- c("Galax urceolata", "Galax aphylla")
Pyxidanthera_barbulata <- c("Pyxidanthera barbulata", "Pyxidanthera barbulata var. barbulata")
Pyxidanthera_brevifolia <- c("Pyxidanthera brevifolia", "Pyxidanthera barbulata var. brevifolia")

# Use gatoRs to download records with resolved synonyms
gators_download(synonyms.list = Shortia_galacifolia,
                write.file = TRUE,
                filename = "data/01_download/raw/Shortia_galacifolia_raw_2025_06_27.csv")

gators_download(synonyms.list = Galax_urceolata,
                write.file = TRUE,
                filename = "data/01_download/raw/Galax_urceolata_raw_2025_06_27.csv")

gators_download(synonyms.list = Pyxidanthera_barbulata,
                write.file = TRUE,
                filename = "data/01_download/raw/Pyxidanthera_barbulata_raw_2025_06_27.csv")

gators_download(synonyms.list = Pyxidanthera_brevifolia,
                write.file = TRUE,
                filename = "data/01_download/raw/Pyxidanthera_brevifolia_raw_2025_06_27.csv")

## ---- C) Preview Downloaded Files ----

# Read one downloaded file
rawdf <- read.csv("data/01_download/raw/Shortia_galacifolia_raw_2025_06_27.csv")

# Preview columns and dimensions
names(rawdf)
nrow(rawdf)

# Visualize records interactively - The error message here indicates many points do not have long/lat values 
leaflet(rawdf) %>% 
  addTiles() %>% 
  addMarkers(label = paste0(rawdf$longitude, ", ", rawdf$latitude))
