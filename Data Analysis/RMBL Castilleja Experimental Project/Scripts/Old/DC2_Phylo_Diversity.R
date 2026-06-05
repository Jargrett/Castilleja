install.packages("V.PhyloMaker2")
library(V.PhyloMaker2)

# Create a species list in the required format
# Needs columns: species, genus, family
species_list <- data.frame(
  species = c("Fragaria_virginiana", "Erigeron_coulteri"), # all your species
  genus   = c("Fragaria", "Erigeron"),
  family  = c("Rosaceae", "Asteraceae")
)

# Generate phylogeny
phylo <- phylo.maker(sp.list = species_list, 
                     tree    = GBOTB.extended.TPL,
                     nodes   = nodes.info.1,
                     output.tree = TRUE)

# Extract the phylogenetic tree
tree <- phylo$scenario.3  # scenario 3 is generally recommended

install.packages("picante")
library(picante)

# Your species matrix with plots as rows and species as columns
# Species names must match the tree tip labels exactly
species_matrix <- presence.data %>% select(17:158)

# Match species matrix to tree
pruned <- match.phylo.comm(tree, species_matrix)

# Faith's Phylogenetic Diversity (PD) — most common metric
pd_results <- pd(pruned$comm, pruned$phy, include.root = TRUE)

# Mean Pairwise Distance (MPD) — measures average relatedness
mpd_results <- mpd(pruned$comm, cophenetic(pruned$phy))

# Mean Nearest Taxon Distance (MNTD) — measures relatedness of closest relatives
mntd_results <- mntd(pruned$comm, cophenetic(pruned$phy))

# Add to your dataframe
presence.data$PD   <- pd_results$PD
presence.data$MPD  <- mpd_results
presence.data$MNTD <- mntd_results

library(lme4)
library(lmerTest)

# Faith's PD
mod_pd <- lmer(PD ~ castilleja + (1|pair),
               data = castilleja.cover)  # use full dataset not just presence
summary(mod_pd)

# MPD
mod_mpd <- lmer(MPD ~ castilleja + (1|pair),
                data = castilleja.cover)
summary(mod_mpd)

ggplot(castilleja.cover, aes(x = castilleja, y = PD, fill = castilleja)) +
  geom_boxplot(alpha = 0.7) +
  labs(x = "Castilleja Presence",
       y = "Faith's Phylogenetic Diversity") +
  theme_classic() +
  theme(legend.position = "none")