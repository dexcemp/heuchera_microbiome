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

### mixed model approach -- Faith PD -- favors null model with no predictors
# manual modeling approach
library(lme4)
model <- lmer(faith_pd ~ (1 | host_species) + latitude + BIOCLIM_1 + BIOCLIM_12 + aridity_index_UNEP + GTOPO30_ELEVATION + phX10 + sand + carbon, data = combined.scaled.reduced)
model.norandom <- lm(faith_pd ~ latitude + BIOCLIM_1 + BIOCLIM_12 + aridity_index_UNEP + GTOPO30_ELEVATION + phX10 + sand + carbon, data = combined.scaled.reduced)
model.latitude <- lm(faith_pd ~ latitude, data = combined.scaled.reduced)
model.BIO1 <- lm(faith_pd ~ BIOCLIM_1, data = combined.scaled.reduced)
model.BIO12 <- lm(faith_pd ~ BIOCLIM_12, data = combined.scaled.reduced)
model.phX10 <- lm(faith_pd ~ phX10, data = combined.scaled.reduced)
model.aridity <- lm(faith_pd ~ aridity_index_UNEP, data = combined.scaled.reduced)
model.simple.random <- lmer(faith_pd ~ (1 | host_species), data = combined.scaled.reduced)
model.null <- lm(faith_pd ~ 1, data = combined.scaled.reduced)
# Calculate Akaike weights
library(qpcR)
format(akaike.weights(c(AIC(model), AIC(model.norandom), AIC(model.latitude), AIC(model.BIO1), AIC(model.BIO1), AIC(model.BIO12), AIC(model.phX10), AIC(model.aridity), AIC(model.simple.random), AIC(model.null))), scientific=F)

# automated modeling approach
model.full.faith_pd <- lmer(faith_pd ~ (1 | host_species) + latitude + BIOCLIM_1 + BIOCLIM_12 + aridity_index_UNEP + GTOPO30_ELEVATION + phX10 + sand + carbon, data = combined.scaled.reduced)
library(lmerTest)
step(model.full.faith_pd)

### mixed model approach -- Shannon Diversity -- favors null model with no predictors
model.full.shannon <- lmer(shannon_entropy ~ (1 | host_species) + latitude + BIOCLIM_1 + BIOCLIM_12 + aridity_index_UNEP + GTOPO30_ELEVATION + phX10 + sand + carbon, data = combined.scaled.reduced)
step(model.full.shannon)


### distance matrix tests
# Geographic distance matrix on full sample set
library(geodist)
geodistances <- geodist_vec(combined$longitude, combined$latitude, measure = "geodesic")
# Import QIIME distance matrix 
unifracdistance <- read.table("distance-matrix.tsv", header = TRUE, row.names = 1, check.names=FALSE)
# Partial Mantel test
library(ape)
mantel.test(geodistances, unifracdistance, nperm = 10000)

# Load tree of all samples and nothing but the samples
tree <- read.tree("NuclearTreePartitionedByGene.rooted.FINAL.tre")
treedist <- cophenetic(tree)
# Make sure ordering matches the qiime distance matrix
treedist <- treedist[order(rownames(unifracdistance)), order(colnames(unifracdistance))] 
mantel.test(treedist, unifracdistance, nperm = 10000)

soildist <- dist(combined[,c("phX10", "sand", "carbon")], upper = TRUE)
soildist <- as.data.frame(as.matrix(soildist))
mantel.test(soildist, unifracdistance, nperm = 10000)
## Partial mantel test if needed
#library(vegan)
#mantel.partial(soildist, unifracdistance, geodistances, permutations = 10000)

environmentdist <- dist(combined[,c("latitude", "BIOCLIM_1", "BIOCLIM_12", "GTOPO30_ELEVATION")], upper = TRUE)
environmentdist <- as.data.frame(as.matrix(environmentdist))
mantel.test(environmentdist, unifracdistance, nperm = 10000)
## Partial mantel test if needed
#mantel.partial(environmentdist, unifracdistance, geodistances, permutations = 10000)


### phylogenetic signal tests... asking if diversity is related to plant phylogeny (close relatives maintain similar diversity)
library(phytools)
faith_pd_vector <- combined$faith_pd
names(faith_pd_vector) <- combined$sample.id
phylosig(tree, faith_pd_vector, method = "lambda", test = TRUE)
shannon_vector <- combined$shannon_entropy
names(shannon_vector) <- combined$sample.id
phylosig(tree, shannon_vector, method = "lambda", test = TRUE)
