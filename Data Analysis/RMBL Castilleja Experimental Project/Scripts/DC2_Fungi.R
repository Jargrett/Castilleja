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

devtools::install_github("brendanf/FUNGuildR")
library(lme4)
library(lmerTest)  # provides p-values for lmer objects

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
otu_key <- data.frame(raw_id = otu_tax$OTU_ID, OTU_ID  = paste0("OTU", seq_len(nrow(otu_tax))))

# TAXONOMY QA/QC
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
tax_clean <- Taxonomy %>% filter(!is.na(Phylum) & Phylum != "")#from 4111 to 3875

droped.tax   <- which(is.na(Taxonomy$Phylum) | Taxonomy$Phylum == "")#saved which taxa were dropped


# Rename hash to OTU1, 2, 3, 4, ect...
otu_tax_id <- otu_tax %>%
  rename(raw_id = OTU_ID) %>%
  left_join(otu_key, by = "raw_id") %>%
  select(OTU_ID, everything(), -raw_id, -confidence, -taxonomy)

# Apply to taxonomy table
Taxonomy <- Taxonomy %>%
  rename(raw_id = OTU_ID) %>%
  left_join(otu_key, by = "raw_id")

tax_clean <- tax_clean %>%
  rename(raw_id = OTU_ID) %>%
  left_join(otu_key, by = "raw_id")

# Get the OTU IDs that survived the taxonomy filter in Taxonomy1
kept_otu <- tax_clean$OTU_ID

# Filter otu_tax_id to only include those OTUs
otu_tax_id <- otu_tax_id %>% filter(OTU_ID %in% kept_otu)

#create OTU matrix
otu_matrix <- otu_tax_id %>%
  column_to_rownames("OTU_ID") %>%  # set OTU_ID as row names
  t() %>% # transpose: samples become rows, OTUs become columns
  as.data.frame()  

#Remove Sample EL37A
otu_matrix <- otu_matrix[otu_matrix$OTU1 != 90, ]

metadata <- metadata %>% 
  dplyr::rename("sample_id" = "sample-id") %>% 
  filter (sample_id != "EL37A")

# RAREFACTION 
# ---------------------------------------------------------------------------
raremax <- min(rowSums(otu_matrix))#Min depth is at 49490 <- this means that the min was not set by EL37A
max(rowSums(otu_matrix))#Max depth is at 168870
#Rarefaction curves based on raw sequence data--this takes forever to run
rarecurve(otu_matrix, step = 100, sample = raremax, label=F, 
          main = "Pre-rarefaction") #Mostly looks pretty reasonable and saturating
small.samples <- otu_matrix[,which(rowSums(otu_matrix)<150000)]
rarecurve(small.samples,step=100,label=T, main = "Small Samples") #Look specifically at samples that have <150K depth
#Double check this intuition by comparing observed OTUs to rarefied richness
metadata$Obs.Rich <- specnumber(otu_matrix)
metadata$Rare.Rich <- rarefy(otu_matrix, raremax)
ggscatter(data = metadata, 'Obs.Rich', 'Rare.Rich', color = 'litter', shape = 'removal',
          xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
#pretty much on 1:1 line

#However, sequencing depth does drive observed diversity patterns
metadata$Seqs<-rowSums(otu_matrix)
ggscatter(metadata,'Seqs','Obs.Rich', add = "reg.line",color = 'litter', facet.by = c('removal'))
#The relationship between sequence depth and richness depends on treatment, and is mostly positive
ggscatter(metadata,'Seqs','Rare.Rich', add = "reg.line", color = 'litter', facet.by = c('removal'))

#Rarefy to min seq depth for richness analyses
metadata$Rare.Rich <- rarefy(otu_matrix, raremax)
rarefied <- rrarefy(otu_matrix, raremax)
otu_rarefied <- as.data.frame(rarefied)

otu_matrix_names_rare <- tibble::rownames_to_column(otu_rarefied, "sample_id")

#COMPOSITION
# ---------------------------------------------------------------------------
set.seed(123)

mean.avg.dist <- as.matrix(avgdist(otu_matrix, sample = raremax, iterations = 100))
summary(mean.avg.dist)
hist(mean.avg.dist)

fungi.nmds <- metaMDS(mean.avg.dist, k = 3, try = 500)
fungi.nmds # stress at k = 3: 0.15
stressplot(fungi.nmds)

nmds.scores <- as.data.frame(scores(fungi.nmds))
nmds.scores <- tibble::rownames_to_column(nmds.scores, "sample_id")

NMDS <- metadata %>% 
  left_join(nmds.scores, by = "sample_id")

# NMDS plot -- colored by litter, shaped by removal
y1.NMDS <- filter(NMDS, year == "2024")

y1.fungal <- ggplot(y1.NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = litter, shape = removal), size = 3) +
  theme_classic() +
  theme(legend.position = "right", legend.text = element_text(size = 12),
        legend.title = element_text(size = 14), axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  coord_fixed() +
  labs(colour = "Litter", shape = "Removal",
       title  = "Fungal Composition (2024)")

y2.NMDS <- filter(NMDS, year == "2025")

y2.fungal <- ggplot(y2.NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = litter, shape = removal), size = 3) +
  theme_classic() +
  theme(legend.position = "right", legend.text = element_text(size = 12),
        legend.title = element_text(size = 14), axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  coord_fixed() +
  labs(colour = "Litter", shape = "Removal",
       title  = "Fungal Composition (2025)")


nmds.fungi <- ggarrange(y1.fungal, y2.fungal,
                         nrow = 1, common.legend = T,legend = "right")
nmds.fungi
#Multivariate
NMDS <- NMDS[match(rownames(mean.avg.dist), NMDS$sample_id), ]
  
cap.mod <- capscale(mean.avg.dist ~ removal*litter*year + Condition(block), 
                    data = NMDS, )

perms <- how(blocks = NMDS$block, nperm = 9999)

anova_res <- anova(cap.mod, permutations = perms, by = "terms")

#---------------Diversity Calculations---------------#
library(vegan)

# Average alpha diversity across 100 rarefaction iterations
set.seed(123)
n_iter <- 100

# Initialize storage
shannon_matrix  <- matrix(NA, nrow = nrow(otu_matrix), ncol = n_iter)
richness_matrix <- matrix(NA, nrow = nrow(otu_matrix), ncol = n_iter)
evenness_matrix   <- matrix(NA, nrow = nrow(otu_matrix), ncol = n_iter)

for (i in 1:n_iter) {
  rare_i <- rrarefy(otu_matrix, raremax)
  shannon_matrix[, i]  <- diversity(rare_i, index = "shannon")
  richness_matrix[, i] <- specnumber(rare_i)
  evenness_matrix[, i]   <- diversity(rare_i, index = "shannon") / log(specnumber(rare_i))
}

# Average across iterations
alpha_df <- data.frame(
  sample_id = rownames(otu_matrix),
  shan_rare   = rowMeans(shannon_matrix),
  rich_rare  = rowMeans(richness_matrix),
  even_rare  = rowMeans(evenness_matrix)
) %>%
  left_join(metadata, by = "sample_id")


#shannon
div.lmm <- lmer(shan_rare ~ removal*litter*year + (1|block) + (1|pair), data = alpha_df)
summary(div.lmm)
Anova(div.lmm)# year p = 0.0001783, removal:litter p = 0.0528
emmeans(div.lmm, pairwise ~ removal|litter)
emmip(div.lmm, litter ~ removal)

#richness
rich.lmm <- lmer(rich_rare ~ removal*litter*year + (1|block) + (1|pair), data = alpha_df)
summary(rich.lmm)
Anova(rich.lmm)# year p = 5.827e-10, removal:litter p = 0.0045
emmeans(rich.lmm, pairwise ~ removal|litter)
emmip(rich.lmm, litter ~ removal)

#evenness
even.lmm <- lmer(even_rare ~ removal*litter*year + (1|block) + (1|pair), data = alpha_df)
summary(even.lmm)
Anova(even.lmm)# year = 0.0422
emmeans(even.lmm, pairwise ~ removal|year)
emmip(even.lmm, removal ~ year)


div.mean <- alpha_df %>% 
  group_by(year, removal, litter) %>% 
  dplyr::summarise(mean = mean(shan_rare),
                   se = sd(shan_rare)/sqrt(n()))

rich.mean <- alpha_df %>% 
  group_by(year, removal, litter) %>% 
  dplyr::summarise(mean = mean(rich_rare),
                   se = sd(rich_rare)/sqrt(n()))

even.mean <- alpha_df %>% 
  group_by(year, removal, litter) %>% 
  dplyr::summarise(mean = mean(even_rare),
                   se = sd(even_rare)/sqrt(n()))


fun.rich.plot <- ggplot(data = rich.mean, aes(x = removal, y = mean)) +
  geom_point(data = alpha_df, aes(x = removal, y = rich_rare, color = litter),
             position = position_jitterdodge(0.3, dodge.width = .3), size = 2, alpha = 0.5) +
  geom_point(aes(shape = removal, color = litter), shape = 18, size = 4.5, 
             position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, color = litter),
                position = position_dodge(width = 0.3), width = 0.07) +
  geom_line(aes(color = litter), 
            position = position_dodge(width = 0.3)) +
  theme_pubr() +
  facet_wrap(~year) +
  theme(strip.text = element_text(size = 15),
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent",
                                    color = "gray", linewidth = 0.12)) +
  ylim(0, 555) +
  labs(x = "Castilleja", y = "OTU Richness")

fun.rich.plot

