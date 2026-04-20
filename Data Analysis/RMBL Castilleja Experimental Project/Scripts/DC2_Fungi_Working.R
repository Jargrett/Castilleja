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
library(lmerTest)  # provides p-values for lmer objects

# ---------------------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------------------
# OTU table merged with metadata (one row per sample)
otu_meta <- read_csv("Raw Data/otu_table_with_metadata.csv")

# OTU table with taxonomy
otu_tax <- read_csv("Raw Data/otu_table_with_taxonomy.csv")

# Metadata only
metadata <- read_csv("Raw Data/metadata.csv")

# Ensure year and treatment columns are factors
metadata$year    <- as.factor(metadata$year)
metadata$litter  <- as.factor(metadata$litter)
metadata$removal <- as.factor(metadata$removal)
metadata$block   <- as.factor(metadata$block)
metadata$pair    <- as.factor(metadata$pair)

# Taxonomy only -- for joining later
Taxonomy <- otu_tax %>%
  select(OTU_ID, taxonomy, confidence) %>%
  # Split UNITE taxonomy string into ranks
  # UNITE format: k__Fungi;p__Basidiomycota;c__...
  separate(taxonomy,
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep  = ";",
           fill = "right") %>%
  mutate(across(Kingdom:Species, ~ str_trim(str_remove(., "^[kpcofgs]__"))))

# Filter out OTUs with no phylum assignment
Taxonomy1 <- Taxonomy %>% filter(!is.na(Phylum) & Phylum != "")
dropTAX   <- which(is.na(Taxonomy$Phylum) | Taxonomy$Phylum == "")
# ---------------------------------------------------------------------------
# 1b. RENAME OTU IDs TO OTU1, OTU2, OTU3...
# ---------------------------------------------------------------------------
# QIIME2 uses hashed sequence IDs by default -- rename to readable OTU labels.
# Original hash IDs are kept in hash_id column for reference.

# Build a lookup table: hash_id -> OTU_ID
otu_key <- data.frame(
  hash_id = otu_tax$OTU_ID,
  OTU_ID  = paste0("OTU", seq_len(nrow(otu_tax)))
)

# Apply to taxonomy table
Taxonomy <- Taxonomy %>%
  rename(hash_id = OTU_ID) %>%
  left_join(otu_key, by = "hash_id")

Taxonomy1 <- Taxonomy1 %>%
  rename(hash_id = OTU_ID) %>%
  left_join(otu_key, by = "hash_id")

# Apply to otu_meta -- rename OTU columns from hash to OTU1, OTU2...
hash_cols <- intersect(names(metadata), otu_key$hash_id)
rename_vec <- setNames(hash_cols, otu_key$OTU_ID[match(hash_cols, otu_key$hash_id)])
metadata <- metadata %>% rename(any_of(rename_vec))

# Apply to OTUrarefied after rarefaction (done below in Step 3)
# -- handled automatically since otu_meta column names are already updated




# ---------------------------------------------------------------------------
# 2. PREPARE OTU MATRIX
# ---------------------------------------------------------------------------
# Separate metadata columns from OTU count columns
meta_cols <- c("sample-id", "year", "year_code", "site", "site_code",
               "plot", "pair", "block", "litter", "removal", "removal_held",
               "elevation", "disturbance_distance", "disturbance_type",
               "soil_moisture", "summer_precip", "summer_temp",
               "snowpack", "snowfall", "plate_row", "plate_column",
               "project_id", "ggbc_filename")

otu_meta <- as.data.frame(otu_meta)
OTUonly <- otu_meta %>%
  select(-any_of(meta_cols)) %>%
  select(where(is.numeric))


rownames(OTUonly) <- otu_meta$`sample-id`

# Drop OTUs with zero total abundance
# Drop OTUs with zero total abundance
dropOTU <- which(colSums(OTUonly) == 0)
cat("Zero-abundance OTUs to drop:", length(dropOTU), "\n")

if (length(dropOTU) > 0) {
  OTUonly1 <- OTUonly[, -dropOTU]
} else {
  OTUonly1 <- OTUonly
}

cat("OTUs remaining:", ncol(OTUonly1), "\n")
# Drop OTUs with no taxonomy
# Drop OTUs with no taxonomy
valid_dropTAX <- dropTAX[dropTAX <= ncol(OTUonly1)]

if (length(valid_dropTAX) > 0) {
  OTUonly2 <- OTUonly1[, -valid_dropTAX]
} else {
  OTUonly2 <- OTUonly1
}

cat("OTUs after taxonomy filter:", ncol(OTUonly2), "\n")

# ---------------------------------------------------------------------------
# 3. RAREFY
# ---------------------------------------------------------------------------
mat     <- as.matrix(OTUonly2)
raremax <- min(rowSums(mat))
cat("Rarefying to:", raremax, "reads\n")

# Rarefaction curves -- check that curves plateau before raremax
rarecurve(mat, step = 100, sample = raremax, col = "blue", label = FALSE,
          main = "Rarefaction curves -- pre-rarefaction")

rarefied    <- rrarefy(mat, raremax)
rarecurve(rarefied, step = 100, sample = raremax, col = "blue", label = FALSE,
          main = "Rarefaction curves -- post-rarefaction")

OTUrarefied <- as.data.frame(rarefied)
OTUrarefied <- rownames_to_column(OTUrarefied, "sample-id")
# OTU columns are now named OTU1, OTU2, OTU3... matching otu_key lookup table

# Identify low-diversity samples after rarefaction
low_samples <- rownames(OTUrarefied)[rowSums(OTUrarefied[,-1]) < 10000]
print(low_samples)

# Check their metadata
metadata %>% filter(`sample-id` %in% low_samples)

# ---------------------------------------------------------------------------
# 4. BETA DIVERSITY -- NMDS + PERMANOVA (all samples)
# ---------------------------------------------------------------------------
set.seed(123)

EL.dist <- vegdist(OTUrarefied[, -1], method = "bray")

EL.nmds <- metaMDS(EL.dist, distance = "bray", k = 3, try = 500, trymax = 100)
EL.nmds
stressplot(EL.nmds)

nmds.scores <- as.data.frame(scores(EL.nmds))
NMDS     <- cbind(metadata, nmds.scores)

# NMDS plot -- colored by litter, shaped by removal
ggplot(EL.info, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = litter, shape = removal), size = 3) +
  theme_classic() +
  theme(legend.position    = "right",
        legend.text        = element_text(size = 12),
        legend.title       = element_text(size = 14),
        axis.text          = element_text(size = 14),
        axis.title         = element_text(size = 16)) +
  coord_fixed() +
  labs(colour = "Litter", shape = "Removal",
       title  = "Fungal community composition")



#Multivariate
cap.mod <- capscale(EL.dist ~ removal*litter*year + disturbance_distance + Condition(pair), 
                    data = NMDS)

perms <- how(blocks = NMDS$pair, nperm = 9999)

anova_res <- anova(cap.mod, permutations = perms, by = "terms")


# ---------------------------------------------------------------------------
# 5. BETA DIVERSITY BY SUBSET -- litter treatment
# ---------------------------------------------------------------------------
# Repeat NMDS and PERMANOVA within each litter treatment group
# to test whether removal effect holds within each litter type

for (litter_type in c("Mixed", "Castilleja", "Community", "Control")) {
  
  idx     <- which(EL.info$litter == litter_type)
  sub_otu <- OTUrarefied[idx, -1]
  sub_meta <- EL.info[idx, ]
  
  sub.dist <- vegdist(sub_otu, method = "bray")
  cat("\n--- Litter:", litter_type, "---\n")
  print(adonis2(sub.dist ~ removal, data = sub_meta, permutations = 999))
}


# ---------------------------------------------------------------------------
# 6. ALPHA DIVERSITY
# ---------------------------------------------------------------------------
OTUmat <- as.matrix(OTUrarefied[, -1])

shannon  <- diversity(OTUmat, index = "shannon")
simpson  <- diversity(OTUmat, index = "simpson")
richness <- specnumber(OTUmat)
H        <- diversity(OTUmat)
pielou_J <- H / log(richness)

alpha_df <- data.frame(
  `sample-id` = OTUrarefied$`sample-id`,
  shannon      = shannon,
  simpson      = simpson,
  richness     = richness,
  pielou_J     = pielou_J
) %>%
  left_join(metadata, by = "sample-id")

# ---------------------------------------------------------------------------
# MIXED MODELS -- alpha diversity
# Fixed:  litter * removal * year (full factorial with all interactions)
# Random: (1|block/pair) -- pair nested within block
# ---------------------------------------------------------------------------
# Ensure year is a factor (2 levels: 2024, 2025)
alpha_df$year    <- as.factor(alpha_df$year)
alpha_df$litter  <- as.factor(alpha_df$litter)
alpha_df$removal <- as.factor(alpha_df$removal)
alpha_df$block   <- as.factor(alpha_df$block)
alpha_df$pair    <- as.factor(alpha_df$pair)

# Shannon diversity

# Richness

# Pielou's J

# Simpson
# ---------------------------------------------------------------------------
# If the three-way interaction is not significant, simplify to:
# litter + removal + year + litter:removal + litter:year + removal:year
# ---------------------------------------------------------------------------
# lmer_shannon_reduced <- lmer(shannon ~ litter + removal + year +
#                                litter:removal + litter:year + removal:year +
#                                (1|block/pair), data = alpha_df)
# anova(lmer_shannon, lmer_shannon_reduced)  # likelihood ratio test

# Boxplots for visualization
boxplot(shannon ~ litter,  data = alpha_df, ylab = "Shannon diversity",
        main = "Shannon by litter treatment")
boxplot(shannon ~ removal, data = alpha_df, ylab = "Shannon diversity",
        main = "Shannon by Castilleja removal")
boxplot(shannon ~ year,    data = alpha_df, ylab = "Shannon diversity",
        main = "Shannon by year")


# ---------------------------------------------------------------------------
# 7. RELATIVE ABUNDANCE
# ---------------------------------------------------------------------------
# Long format: join OTU counts with metadata and taxonomy
OTU_long <- OTUrarefied %>%
  pivot_longer(!`sample-id`, names_to = "OTU_ID", values_to = "Abundance") %>%
  inner_join(metadata,   by = "sample-id") %>%
  inner_join(Taxonomy1,  by = "OTU_ID") %>%
  group_by(`sample-id`) %>%
  mutate(relative_abundance = Abundance / sum(Abundance) * 100) %>%
  ungroup()

# By Phylum -- grouped by litter and removal
rel_abund_phylum <- OTU_long %>%
  group_by(litter, removal, Phylum) %>%
  summarize(mean_rel_abund = mean(relative_abundance), .groups = "drop") %>%
  filter(Phylum != "Unassigned" & !is.na(Phylum))

cols <- c("#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#CC6600",
          "#999999", "#000000", "#E69F00", "#56B4E9", "#009E73")

rel_abund_phylum %>%
  ggplot(aes(x = removal,
             y = mean_rel_abund,
             fill = fct_reorder(Phylum, mean_rel_abund))) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +
  labs(x = "Castilleja removal", y = "Mean Relative Abundance (%)") +
  theme_classic() +
  theme(legend.position  = "right",
        legend.title     = element_blank(),
        axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.text      = element_text(size = 8)) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE)) +
  scale_fill_manual(values = cols)


# ---------------------------------------------------------------------------
# 8. INDICATOR SPECIES
# ---------------------------------------------------------------------------
# By litter treatment
litter_groups <- metadata$litter
inv_litter     <- multipatt(OTUonly2, litter_groups,
                            func = "r.g", control = how(nperm = 999))
summary(inv_litter)

# By removal
removal_groups <- metadata$removal
inv_removal     <- multipatt(OTUonly2, removal_groups,
                             func = "r.g", control = how(nperm = 999))
summary(inv_removal)


# ---------------------------------------------------------------------------
# 9. FUNGUILD TROPHIC MODE ANALYSIS
# ---------------------------------------------------------------------------
assigned <- funguild_assign(OTU_long)
assigned$trophicMode[is.na(assigned$trophicMode)] <- "Unassigned"

confidence_levels <- c("Highly Probable", "Probable")
assigned1 <- assigned %>% filter(confidenceRanking %in% confidence_levels)

# Trophic mode by litter and removal
assigned1 %>%
  ggplot(aes(x = removal,
             y = relative_abundance,
             fill = fct_reorder(trophicMode, relative_abundance))) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +
  facet_grid(~ litter) +
  labs(x = "Castilleja removal", y = "Relative Abundance") +
  theme_classic() +
  theme(legend.position = "right",
        legend.title    = element_blank(),
        axis.text.x     = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.text     = element_text(size = 8)) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE))

# LMs for each trophic mode
guilds <- assigned1 %>%
  select(`sample-id`, litter, removal, block, OTU_ID, relative_abundance, trophicMode)

guilds_wide <- guilds %>%
  group_by(`sample-id`, litter, removal, block, trophicMode) %>%
  summarize(total_abund = sum(relative_abundance), .groups = "drop") %>%
  pivot_wider(names_from = trophicMode, values_from = total_abund, values_fill = 0)

# Test each trophic mode against litter and removal
guilds_wide$year    <- as.factor(guilds_wide$year)
guilds_wide$litter  <- as.factor(guilds_wide$litter)
guilds_wide$removal <- as.factor(guilds_wide$removal)
guilds_wide$block   <- as.factor(guilds_wide$block)
guilds_wide$pair    <- as.factor(guilds_wide$pair)

for (mode in c("Pathotroph", "Symbiotroph", "Saprotroph")) {
  if (mode %in% names(guilds_wide)) {
    cat("\n---", mode, "---\n")
    lmer_mode <- lmer(guilds_wide[[mode]] ~ litter * removal * year +
                        (1|block/pair), data = guilds_wide)
    print(anova(lmer_mode))
  }
}
