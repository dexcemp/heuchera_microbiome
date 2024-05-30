# import data
metadata <- read.table("sample-metadata.tsv", sep = "\t", header = TRUE)
pd <- read.table("faith_pd_vector.tsv", sep = "\t", header = TRUE)
shannon <- read.table("shannon_vector.tsv", sep = "\t", header = TRUE)
combined <- merge(metadata, pd, by = "sample.id")
combined <- merge(combined, shannon, by = "sample.id")

# scale means and variance
library(dplyr)
combined.scaled <- combined %>% mutate_if(is.numeric,scale)

# simplify column names
names(combined.scaled) <- c("sample.id", "host_species", "host_subsection", "latitude", "longitude", "BIOCLIM_1", "BIOCLIM_12", "aridity_index_UNEP", "GTOPO30_ELEVATION", "phX10", "sand", "carbon", "faith_pd", "shannon_entropy")
names(combined) <- c("sample.id", "host_species", "host_subsection", "latitude", "longitude", "BIOCLIM_1", "BIOCLIM_12", "aridity_index_UNEP", "GTOPO30_ELEVATION", "phX10", "sand", "carbon", "faith_pd", "shannon_entropy")

# subset to just the parviflora group
combined.scaled.parviflora <- combined.scaled %>% filter(host_species %in% c("Heuchera_missouriensis", "Heuchera_puberula", "Heuchera_parviflora_parviflora", "Heuchera_parviflora_saurensis"))
combined.parviflora <- combined %>% filter(host_species %in% c("Heuchera_puberula", "Heuchera_parviflora_parviflora", "Heuchera_parviflora_saurensis"))

# ANOVA for parviflora group diversity statistics
summary(aov(shannon_entropy ~ host_species, data = combined.scaled.parviflora))
summary(aov(faith_pd ~ host_species, data = combined.scaled.parviflora))

# subset to the well-sampled species
combined.scaled.reduced <- combined.scaled %>% filter(host_species %in% c("Heuchera_americana_americana", "Heuchera_americana_hirsuticaulis", "Heuchera_richardsonii", "Heuchera_longiflora_longiflora", "Heuchera_longiflora_aceroides", "Heuchera_missouriensis", "Heuchera_puberula", "Heuchera_parviflora_parviflora", "Heuchera_parviflora_saurensis"))
combined.reduced <- combined %>% filter(host_species %in% c("Heuchera_americana_americana", "Heuchera_americana_hirsuticaulis", "Heuchera_richardsonii", "Heuchera_longiflora_longiflora", "Heuchera_longiflora_aceroides", "Heuchera_missouriensis", "Heuchera_puberula", "Heuchera_parviflora_parviflora", "Heuchera_parviflora_saurensis"))

# boxplot figs
# Enforce boxplot label order
combined.reduced$host_species <- factor(combined.reduced $host_species, levels = c('Heuchera_americana_americana','Heuchera_americana_hirsuticaulis', "Heuchera_richardsonii", "Heuchera_longiflora_longiflora", "Heuchera_longiflora_aceroides", "Heuchera_missouriensis", "Heuchera_parviflora_parviflora", "Heuchera_parviflora_saurensis", "Heuchera_puberula"),ordered = TRUE)
library(ggplot2)
library(gridExtra)
# Fix titles for ITS vs 16S
p1 <- ggplot(combined.reduced, aes(factor(host_species), faith_pd)) + stat_boxplot() + coord_flip() + xlab("") + ylab("Faith's PD") + ggtitle("Bacterial Faith's PD vs. taxon")
p2 <- ggplot(combined.reduced, aes(factor(host_species), shannon_entropy)) + stat_boxplot() + coord_flip() + xlab("") + ylab("Shannon diversity") + ggtitle("Bacterial Shannon diversity vs. taxon")
grid.arrange(p1, p2, nrow = 2)

### Faith PD 
# manual modeling approach
library(lme4)
model.lmm <- lmer(faith_pd ~ (1 | host_species) + latitude + BIOCLIM_1 + BIOCLIM_12 + aridity_index_UNEP + GTOPO30_ELEVATION + phX10 + sand + carbon, data = combined.scaled.reduced)
model.norandom <- lm(faith_pd ~ latitude + BIOCLIM_1 + BIOCLIM_12 + aridity_index_UNEP + GTOPO30_ELEVATION + phX10 + sand + carbon, data = combined.scaled.reduced)
model.full <- lm(faith_pd ~ host_species + latitude + BIOCLIM_1 + BIOCLIM_12 + aridity_index_UNEP + GTOPO30_ELEVATION + phX10 + sand + carbon, data = combined.scaled.reduced)
model.latitude <- lm(faith_pd ~ latitude, data = combined.scaled.reduced)
model.BIO1 <- lm(faith_pd ~ BIOCLIM_1, data = combined.scaled.reduced)
model.BIO12 <- lm(faith_pd ~ BIOCLIM_12, data = combined.scaled.reduced)
model.phX10 <- lm(faith_pd ~ phX10, data = combined.scaled.reduced)
model.aridity <- lm(faith_pd ~ aridity_index_UNEP, data = combined.scaled.reduced)
model.hostspecies <- lm(faith_pd ~ host_species, data = combined.scaled.reduced)
model.sand <- lm(faith_pd ~ sand, data = combined.scaled.reduced)
model.carbon <- lm(faith_pd ~ carbon, data = combined.scaled.reduced)
model.elevation <- lm(faith_pd ~ GTOPO30_ELEVATION, data = combined.scaled.reduced)
model.simple.random <- lmer(faith_pd ~ (1 | host_species), data = combined.scaled.reduced)
model.null <- lm(faith_pd ~ 1, data = combined.scaled.reduced)

# Calculate Akaike weights using geiger
library(geiger)
format(aicw(c(AIC(model.lmm), AIC(model.norandom), AIC(model.full), AIC(model.latitude), AIC(model.BIO1), AIC(model.BIO12), AIC(model.phX10), AIC(model.aridity), AIC(model.hostspecies), AIC(model.sand), AIC(model.carbon), AIC(model.elevation), AIC(model.simple.random), AIC(model.null))))

# Automated modeling approach
model.full.faith_pd <- lmer(faith_pd ~ (1 | host_species) + latitude + BIOCLIM_1 + BIOCLIM_12 + aridity_index_UNEP + GTOPO30_ELEVATION + phX10 + sand + carbon, data = combined.scaled.reduced)
library(lmerTest)
step(model.full.faith_pd)

### Shannon Diversity
#manual modeling approach
model.lmm <- lmer(shannon_entropy ~ (1 | host_species) + latitude + BIOCLIM_1 + BIOCLIM_12 + aridity_index_UNEP + GTOPO30_ELEVATION + phX10 + sand + carbon, data = combined.scaled.reduced)
model.norandom <- lm(shannon_entropy ~ latitude + BIOCLIM_1 + BIOCLIM_12 + aridity_index_UNEP + GTOPO30_ELEVATION + phX10 + sand + carbon, data = combined.scaled.reduced)
model.full <- lm(shannon_entropy ~ host_species + latitude + BIOCLIM_1 + BIOCLIM_12 + aridity_index_UNEP + GTOPO30_ELEVATION + phX10 + sand + carbon, data = combined.scaled.reduced)
model.latitude <- lm(shannon_entropy ~ latitude, data = combined.scaled.reduced)
model.BIO1 <- lm(shannon_entropy ~ BIOCLIM_1, data = combined.scaled.reduced)
model.BIO12 <- lm(shannon_entropy ~ BIOCLIM_12, data = combined.scaled.reduced)
model.phX10 <- lm(shannon_entropy ~ phX10, data = combined.scaled.reduced)
model.aridity <- lm(shannon_entropy ~ aridity_index_UNEP, data = combined.scaled.reduced)
model.hostspecies <- lm(shannon_entropy ~ host_species, data = combined.scaled.reduced)
model.sand <- lm(shannon_entropy ~ sand, data = combined.scaled.reduced)
model.carbon <- lm(shannon_entropy ~ carbon, data = combined.scaled.reduced)
model.elevation <- lm(shannon_entropy ~ GTOPO30_ELEVATION, data = combined.scaled.reduced)
model.simple.random <- lmer(shannon_entropy ~ (1 | host_species), data = combined.scaled.reduced)
model.null <- lm(shannon_entropy ~ 1, data = combined.scaled.reduced)

# Calculate Akaike weights using geiger
library(geiger)
format(aicw(c(AIC(model.lmm), AIC(model.norandom), AIC(model.full), AIC(model.latitude), AIC(model.BIO1), AIC(model.BIO12), AIC(model.phX10), AIC(model.aridity), AIC(model.hostspecies), AIC(model.sand), AIC(model.carbon), AIC(model.elevation), AIC(model.simple.random), AIC(model.null))))

# Automated modeling approach
model.full.shannon <- lmer(shannon_entropy ~ (1 | host_species) + latitude + BIOCLIM_1 + BIOCLIM_12 + aridity_index_UNEP + GTOPO30_ELEVATION + phX10 + sand + carbon, data = combined.scaled.reduced)
step(model.full.shannon)


### Distance matrix tests
# Import QIIME distance matrix 
unifracdistance <- read.table("distance-matrix.tsv", header = TRUE, row.names = 1, check.names=FALSE)

# Geographic distance matrix on full sample set
library(geodist)
geodistances <- geodist_vec(combined$longitude, combined$latitude, measure = "geodesic")

#Mantel tests
library(vegan)
mantel(geodistances, unifracdistance, permutations = 10000)

# Load tree of all samples and nothing but the samples
tree <- read.tree("NuclearTreePartitionedByGene.rooted.FINAL.tre")
treedist <- cophenetic(tree)

# Make sure ordering matches the qiime distance matrix
treedist <- treedist[order(rownames(unifracdistance)), order(colnames(unifracdistance))] 
mantel(treedist, unifracdistance, permutations = 10000)

#Soil environment distance
soildist <- dist(combined[,c("phX10", "sand", "carbon")], upper = TRUE)
soildist <- as.data.frame(as.matrix(soildist))
mantel(soildist, unifracdistance, permutations = 10000)

## Partial mantel test if needed
mantel.partial(soildist, unifracdistance, geodistances, permutations = 10000)

#Non-soil environment distance
environmentdist <- dist(combined[,c("latitude", "BIOCLIM_1", "BIOCLIM_12", "GTOPO30_ELEVATION")], upper = TRUE)
environmentdist <- as.data.frame(as.matrix(environmentdist))
mantel(environmentdist, unifracdistance, permutations = 10000)

## Partial mantel test if needed
mantel.partial(environmentdist, unifracdistance, geodistances, permutations = 10000)


### Phylogenetic signal tests... asking if diversity is related to plant phylogeny (close relatives maintain similar diversity)
library(phytools)
faith_pd_vector <- combined$faith_pd
names(faith_pd_vector) <- combined$sample.id
phylosig(tree, faith_pd_vector, method = "lambda", test = TRUE)
shannon_vector <- combined$shannon_entropy
names(shannon_vector) <- combined$sample.id
phylosig(tree, shannon_vector, method = "lambda", test = TRUE)


###Re-subset tree depending on rarefaction cutoff
# Re-subset tree to rarefaction at 50 sequences taxon list
treeraw <- read.tree("NuclearTreePartitionedByGene.rooted.samplelabelsver2.unrarefied.tre")
rarefied50_list <- c("A10-1", "A10-2", "A11-2", "A12-1", "A13-6", "A14", "A15-1", "A15-8", "A21-9", "A26-3", "A30-10", "A34-4", "A35-2", "A36", "A38-1", "A39-2", "A4-5", "A40-2", "A42-3", "A45-1", "A49-8", "A6-1", "A7-1", "A8", "FL-1", "FL-10", "FL-13", "FL-16", "FL-17", "FL-20", "FL-21", "FL-22", "FL-23", "FL-24", "FL-25", "FL-26", "FL-27", "FL-3", "FL-6", "FL-7", "FL-8", "FL-9", "H102-3", "H105-1", "H107-2", "H113-1", "H115-2", "H118-1", "H124-1", "H124-2", "H129-1", "H134-2", "H141-2", "H146-2", "H149-2", "H164-2", "H174-4", "H176-4", "H177-2", "H181-2", "H182-4", "H183-4", "H191-4", "H193-2", "H194-2", "H207-4", "H210-4", "H211-4", "H212-4", "H213-4", "H218-2", "H220-4", "H225-4", "H226-4", "H227-4", "H230-6", "H37-1", "H45-1", "H49", "H50-1", "L12-3", "L18-2", "L23-1", "L23-3")
rarefied50_tree <- keep.tip(treeraw, rarefied50_list)

# new tree for rarefied at 50 distance matrix
write.tree(rarefied50_tree, "NuclearTreePartitionedByGene.rooted.FINALrarefied50.tre")








