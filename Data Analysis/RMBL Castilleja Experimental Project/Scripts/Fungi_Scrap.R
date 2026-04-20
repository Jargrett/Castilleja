setwd("~/Desktop/Castilleja/Data Analysis/RMBL Castilleja Experimental Project")
library(readr)
library(tidyverse)
library(tibble)
library(vegan)
library(ggplot2)
library(ggrepel)
library(car)
library(indicspecies)
library(FUNGuildR)
library(lme4)
library(lmerTest)
library(labdsv)#enables restructuring for ecological analysis


# LOAD DATA
# ---------------------------------------------------------------------------
# OTU table merged with metadata (one row per sample)
otu_meta <- read_csv("Raw Data/otu_table_with_metadata.csv")
# OTU table with taxonomy
otu_tax <- read_csv("Raw Data/otu_table_with_taxonomy.csv")
# Metadata only
metadata <- read_csv("Raw Data/metadata.csv")

metadata$year    <- as.factor(metadata$year)
metadata$litter  <- as.factor(metadata$litter)
metadata$removal <- as.factor(metadata$removal)
metadata$block   <- as.factor(metadata$block)
metadata$pair    <- as.factor(metadata$pair)

# RENAME OTU's
# ---------------------------------------------------------------------------
#create dataframe with new OTU ID and has id only
otu_key <- data.frame(hash_id = otu_tax$OTU_ID, OTU_ID  = paste0("OTU", seq_len(nrow(otu_tax))))

# TAXONOMY CHECK and Clean
# ---------------------------------------------------------------------------
#Splitting taxonomy column 
Taxonomy <- otu_tax %>%
  select(OTU_ID, taxonomy, confidence) %>%
  separate(taxonomy,
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep  = ";",
           fill = "right") %>%
  mutate(across(Kingdom:Species, ~ str_trim(str_remove(., "^[kpcofgs]__"))))

# Filter out OTUs with no phylum assignment
Taxonomy1 <- Taxonomy %>% filter(!is.na(Phylum) & Phylum != "")
dropTAX   <- which(is.na(Taxonomy$Phylum) | Taxonomy$Phylum == "")

# See confidence by deepest assigned rank
Taxonomy %>% mutate(resolved_to = case_when(
  !is.na(Species) & Species != "" ~ "Species",
  !is.na(Genus)   & Genus   != "" ~ "Genus",
  !is.na(Family)  & Family  != "" ~ "Family",
  !is.na(Order)   & Order   != "" ~ "Order",
  !is.na(Class)   & Class   != "" ~ "Class",
  !is.na(Phylum)  & Phylum  != "" ~ "Phylum",
  TRUE ~ "Unassigned")) %>%
  group_by(resolved_to) %>%
  summarize(n = n(), mean_conf  = round(mean(as.numeric(confidence), na.rm = TRUE), 3),
            min_conf   = round(min(as.numeric(confidence),  na.rm = TRUE), 3))

# Build lookup table
otu_key <- data.frame(hash_id = otu_tax$OTU_ID,
                      OTU_ID  = paste0("OTU", seq_len(nrow(otu_tax))))

# Rename hash to OTU1, 2, 3, 4, ect...
otu_tax_id <- otu_tax %>%
  rename(hash_id = OTU_ID) %>%
  left_join(otu_key, by = "hash_id") %>%
  select(OTU_ID, everything(), -hash_id, -confidence, -taxonomy)

# Apply to taxonomy table
Taxonomy <- Taxonomy %>%
  rename(hash_id = OTU_ID) %>%
  left_join(otu_key, by = "hash_id")

tax_clean <- Taxonomy1 %>%
  rename(hash_id = OTU_ID) %>%
  left_join(otu_key, by = "hash_id")

# Get the OTU IDs that survived the taxonomy filter in Taxonomy1
kept_otu <- tax_clean$OTU_ID

# Filter otu_tax_id to only include those OTUs
otu_tax_id <- otu_tax_id %>%
  filter(OTU_ID %in% kept_otu)

#create OTU matrix
otu_matrix <- otu_tax_id %>%
  column_to_rownames("OTU_ID") %>%  # set OTU_ID as row names
  t() %>%                            # transpose: samples become rows, OTUs become columns
  as.data.frame()  


# RAREIFCATION 
# ---------------------------------------------------------------------------
#Rarefaction curves based on raw sequence data--this takes forever to run
rarecurve(otu_matrix, step = 100, sample = raremax, label=F, 
          main = "Pre-rarefaction")) #Mostly looks pretty reasonable and saturating
min(rowSums(otu_matrix))#Min depth is at 49490
max(rowSums(otu_matrix))#Max depth is at 168870
small.samples <- otu_matrix[,which(rowSums(otu_matrix)<150000)]
rarecurve(small.samples,step=10) #Look specifically at samples that have <150K depth
#All of them flatten but EL37A. lets remove it.
otu_matrix <- otu_matrix[otu_matrix$OTU1 != 90, ]
metadata <- metadata %>% 
  dplyr::rename("sample_id" = "sample-id") %>% 
  filter (sample_id != "EL37A")
raremax <- min(rowSums(otu_matrix))#Min depth is at 49490
#Double check this intuition by comparing observed OTUs to rarefied richness
metadata$Obs.Rich <- specnumber(otu_matrix)
metadata$Rare.Rich <- rarefy(otu_matrix, raremax)
ggscatter(data = metadata, 'Obs.Rich', 'Rare.Rich', color = 'litter', shape = 'removal',
          xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
#pretty much on 1:1 line

#However, sequencing depth does drive observed diversity patterns
metadata$Seqs<-rowSums(otu_matrix)
ggscatter(metadata,'Seqs','Obs.Rich',
          add = "reg.line",
          color = 'litter', facet.by = c('removal'))
#The relationship between sequence depth and richness depends on treatment, and is mostly positive
ggscatter(metadata,'Seqs','Rare.Rich',
          add = "reg.line",
          color = 'litter', facet.by = c('removal'))
#The relationships still exist using rarefied richness, but slightly less
#So we should still probably rarefy

#Rarefy to min seq depth for richness analyses
metadata$Rare.Rich <- rarefy(otu_matrix, raremax)
rarefied <- rrarefy(otu_matrix, raremax)
otu_rarefied <- as.data.frame(rarefied)

otu_matrix_names_rare <- tibble::rownames_to_column(otu_rarefied, "sample-id")
# BETA DIVERSITY -- NMDS + PERMANOVA (all samples)
# ---------------------------------------------------------------------------
set.seed(123)


mean.avg.dist <- as.matrix(avgdist(otu_matrix, sample = raremax, iterations = 100))
summary(mean.avg.dist)
hist(mean.avg.dist)

fungi.nmds <- metaMDS(mean.avg.dist, distance = "bray", k = 2, try = 500)
fungi.nmds
stressplot(fungi.nmds)

nmds.scores <- as.data.frame(scores(fungi.nmds))
nmds.scores <- tibble::rownames_to_column(nmds.scores, "sample_id")


NMDS <- cbind(metadata, nmds.scores)

NMDS<- metadata %>% 
  left_join(nmds.scores, by = "sample_id")

write.csv(NMDS,"~/Downloads/NMDS.csv", row.names = FALSE)

# NMDS plot -- colored by litter, shaped by removal
ggplot(NMDS_test, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = block, shape = removal), size = 3) +
  theme_classic() +
  theme(legend.position = "right", legend.text = element_text(size = 12),
        legend.title = element_text(size = 14), axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  stat_ellipse(segments = 20, linetype = 2, alpha = 0.5, aes(group = castilleja)) +
  coord_fixed() +
  labs(colour = "block", shape = "Removal",
       title  = "Fungal community composition")


#Multivariate
cap.mod <- capscale(mean.avg.dist ~ removal*litter*year + Condition(pair), 
                    data = NMDS)

perms <- how(blocks = NMDS$pair, nperm = 9999)

anova_res <- anova(cap.mod, permutations = perms, by = "terms")
