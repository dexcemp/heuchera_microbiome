
# activate conda virtual environment
conda activate qiime2-2020.8 

# Make a manifest file of this format
######################
sample-id	forward-absolute-filepath	reverse-absolute-filepath
A40	$PWD/test_seqs/A40_P1.fastq.gz	$PWD/test_seqs/A40_P2.fastq.gz
######################


# import sequences
# Qiime must have gzipped input files in a folder
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.txt --output-path import.qza --input-format PairedEndFastqManifestPhred33V2

# Trimming is to remove primers, set to 16S V3/V4 primers 515F/806R
qiime dada2 denoise-paired --i-demultiplexed-seqs import.qza --p-trim-left-f 19 --p-trim-left-r 20 --p-trunc-len-f 250 --p-trunc-len-r 200 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza

# Taxonomic analysis
qiime feature-classifier classify-sklearn --i-classifier ./greengenes_13_8_otus/classifier_gg_16S.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

# Remove chloroplast and mitochondrial 16S
qiime taxa filter-table --p-exclude c__Chloroplast,f__mitochondria --i-table table.qza --i-taxonomy taxonomy.qza --o-filtered-table table.organellesremoved.qza

qiime taxa filter-seqs --p-exclude c__Chloroplast,f__mitochondria --i-sequences rep-seqs.qza --i-taxonomy taxonomy.qza --o-filtered-sequences rep-seqs.organellesremoved.qza


# Taxa barplots, no organelles
qiime taxa barplot --i-table table.organellesremoved.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.tsv --o-visualization taxa-bar-plots.qzv

## Build GreenGenes database for 16S (this was already done on BotBot) -- these steps take a LONG time -- about 2 hours
#conda activate qiime2-2020.8 
## Format greengenes database (this was already done on BotBot)
#qiime tools import --type 'FeatureData[Sequence]' --input-path ./greengenes_13_8_otus/rep_set/97_otus.fasta --output-path 97_otus.qza
## Format greengenes database (this was already done on BotBot)
#qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path ./greengenes_13_8_otus/taxonomy/97_otu_taxonomy.txt --output-path 97_ref-taxonomy.qza
## Extract representative reads for OTUS from greengenes; set to 16S V3/V4 primers 515F/806R
## Already done on BotBot but could change if settings or primers change -- consider redoing but it takes a long time
#qiime feature-classifier extract-reads --i-sequences 97_otus.qza   --p-f-primer GTGCCAGCMGCCGCGGTAA --p-r-primer GGACTACHVGGGTWTCTAAT --p-trunc-len 250 --p-min-length 200 --p-max-length 400 --o-reads 97_ref-seqs.qza
## Bayesian approach to taxonomic classifications on greengenes (explained in QIIME documentation -- simply taking the closest sequence has issues
#qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads 97_ref-seqs.qza --i-reference-taxonomy 97_ref-taxonomy.qza --o-classifier classifier_gg_16S.qza

# Run standard feature tables
qiime feature-table summarize --i-table table.organellesremoved.qza --o-visualization table.qzv --m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs --i-data rep-seqs.organellesremoved.qza --o-visualization rep-seqs.qzv

qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv

qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.organellesremoved.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza


## Generate sub-group tables for species group tests
qiime feature-table filter-samples --i-table table.organellesremoved.qza --m-metadata-file americana_samples.tsv --o-filtered-table table.organellesremoved.americanagroup.qza
qiime feature-table filter-samples --i-table table.organellesremoved.qza --m-metadata-file longiflora_samples.tsv --o-filtered-table table.organellesremoved.longifloragroup.qza
qiime feature-table filter-samples --i-table table.organellesremoved.qza --m-metadata-file parviflora_samples.tsv --o-filtered-table table.organellesremoved.parvifloragroup.qza

# Diversity statistics and significance
# p-sampling-depth: change per experiment; look across your samples for feature counts -- choose a number as high as possible but lower thant the lowest sample
# Will not work for too few samples
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.organellesremoved.qza --p-sampling-depth 50 --output-dir core-metrics-results --m-metadata-file sample-metadata.tsv 
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/faith_pd_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization core-metrics-results/faith-pd-group-significance.qzv
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/evenness_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization core-metrics-results/evenness-group-significance.qzv
# Due to the longer time needed for beta diversity significance, a specific metadata column is specified to reduce computational time
qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column host_species --o-visualization core-metrics-results/unweighted-unifrac-host_species-significance.qzv --p-pairwise
qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column host_subsection --o-visualization core-metrics-results/unweighted-unifrac-host_subsection-significance.qzv --p-pairwise


# Sub-group diversity statistics
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.organellesremoved.americanagroup.qza --p-sampling-depth 50 --output-dir core-metrics-results_americana --m-metadata-file sample-metadata.tsv
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.organellesremoved.longifloragroup.qza --p-sampling-depth 50 --output-dir core-metrics-results_longiflora --m-metadata-file sample-metadata.tsv
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.organellesremoved.parvifloragroup.qza --p-sampling-depth 50 --output-dir core-metrics-results_parviflora --m-metadata-file sample-metadata.tsv
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results_americana/faith_pd_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization core-metrics-results_americana/faith-pd-group-significance.qzv
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results_longiflora/faith_pd_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization core-metrics-results_longiflora/faith-pd-group-significance.qzv
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results_parviflora/faith_pd_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization core-metrics-results_parviflora/faith-pd-group-significance.qzv
qiime diversity beta-group-significance --i-distance-matrix core-metrics-results_americana/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column host_species --o-visualization core-metrics-results_americana/unweighted-unifrac-host_species-significance.qzv --p-pairwise
qiime diversity beta-group-significance --i-distance-matrix core-metrics-results_longiflora/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column host_species --o-visualization core-metrics-results_longiflora/unweighted-unifrac-host_species-significance.qzv --p-pairwise
qiime diversity beta-group-significance --i-distance-matrix core-metrics-results_parviflora/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column host_species --o-visualization core-metrics-results_parviflora/unweighted-unifrac-host_species-significance.qzv --p-pairwise


# Ordination plots
qiime emperor plot --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza --m-metadata-file sample-metadata.tsv --p-custom-axes BIOCLIM_1 --o-visualization core-metrics-results/unweighted-unifrac-emperor-BIO1.qzv
qiime emperor plot --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza --m-metadata-file sample-metadata.tsv --p-custom-axes BIOCLIM_12 --o-visualization core-metrics-results/unweighted-unifrac-emperor-BIO12.qzv
qiime emperor plot --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza --m-metadata-file sample-metadata.tsv --p-custom-axes ISRICSOILGRIDS_new_average_phx10percent_reduced --o-visualization core-metrics-results/unweighted-unifrac-emperor-pH.qzv
#qiime emperor plot --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza --m-metadata-file sample-metadata.tsv --p-custom-axes host_species --o-visualization core-metrics-results/bray-curtis-emperor-host_species.qzv

# Rarefaction plots
# p-sampling-depth: change per experiment; look across your samples for feature counts -- choose a number as high as possible but lower thant the lowest sample
qiime diversity alpha-rarefaction --i-table table.organellesremoved.qza --i-phylogeny rooted-tree.qza --p-max-depth 50 --m-metadata-file sample-metadata.tsv --o-visualization alpha-rarefaction.qzv



# ANCOM -- differential abundance analysis
# Helps find which genera, families, etc. differ across metadata categories
qiime feature-table filter-samples --i-table table.organellesremoved.qza --m-metadata-file sample-metadata.tsv --p-where "[host_species]='Heuchera_americana_americana'" --o-filtered-table host_species-table.qza
# Imputation method for zero-abundance samples
qiime composition add-pseudocount --i-table host_species-table.qza --o-composition-table comp-host_species-table.qza
qiime composition ancom --i-table comp-host_species-table.qza --m-metadata-file sample-metadata.tsv --m-metadata-column host_species --o-visualization ancom-host_species.qzv

# ANCOM with taxonomic level collapsing -- here level 6 genus
qiime taxa collapse --i-table host_species-table.qza --i-taxonomy taxonomy.qza --p-level 6 --o-collapsed-table host_species-table-l6.qza
qiime composition add-pseudocount --i-table host_species-table-l6.qza --o-composition-table comp-host_species-table-l6.qza
qiime composition ancom --i-table comp-host_species-table-l6.qza --m-metadata-file sample-metadata.tsv --m-metadata-column host_species --o-visualization l6-ancom-host_species.qzv

# adonis
qiime diversity adonis --i-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --p-formula "host_species+host_subsection" --p-permutations 10000 --p-n-jobs 4 --o-visualization ./core-metrics-results/adonis.unweighted_unifrac.qzv


qiime taxa collapse --i-table table.organellesremoved.qza --o-collapsed-table genus_collapse.table.qza --p-level 6 --i-taxonomy taxonomy.qza
qiime feature-table relative-frequency --i-table genus_collapse.table.qza --o-relative-frequency-table genus_collapse.frequency.table.qza --output-dir genus_collapse.frequency
qiime tools export --input-path genus_collapse.frequency.table.qza --output-path genus_collapse.frequency
biom convert -i ./genus_collapse.frequency/feature-table.biom -o genus_collapse.frequency.table.txt --header-key “taxonomy” --to-tsv
rm -rf genus_collapse.frequency
# Then manually format by adding top two tsv columns for class and subclass


# Analyse how much was thrown away as host DNA using new taxonomic barplots
qiime feature-classifier classify-sklearn --i-classifier ./greengenes_13_8_otus/classifier_gg_16S.qza --i-reads rep-seqs.qza --o-classification taxonomy.withhost.qza
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.withhost.qzv
qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.withhost.qza --m-metadata-file sample-metadata.tsv --o-visualization taxa-bar-plots.withhost.qzv



#######
# Aggregate taxon barplots to host species

qiime feature-table group --i-table table.hostremoved.qza --m-metadata-file sample-metadata.tsv --m-metadata-column host_species --p-axis sample --p-mode mean-ceiling --o-grouped-table table.grouped_hostspecies # note lack of extension
qiime taxa barplot --i-table table.grouped_hostspecies.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.grouped_hostspecies.tsv --o-visualization taxa-bar-plots.grouped_hostspecies.qzv




#######
# Export tsv

qiime tools export --input-path core-metrics-results/faith_pd_vector.qza --output-path export
mv export/alpha-diversity.tsv export/faith_pd_vector.tsv
qiime tools export --input-path core-metrics-results/shannon_vector.qza --output-path export
mv export/alpha-diversity.tsv export/shannon_vector.tsv
# Add sample column header
sed -i 's/^\t/sample-id\t/g' export/*.tsv

qiime tools export --input-path core-metrics-results/unweighted_unifrac_distance_matrix.qza --output-path export
mv export/alpha-diversity.tsv export/unweighted_unifrac_distance_matrix.tsv



#######
# PICRUST
conda activate qiime2-2020.8
qiime tools export --input-path rep-seqs.qza --output-path export
mv ./export/dna-sequences.fasta ./export/rep-seqs.fasta
qiime tools export --input-path table.qza --output-path export
mv ./export/feature-table.biom ./export/table.biom

conda activate picrust2
picrust2_pipeline.py -s ./export/rep-seqs.fasta -i ./export/table.biom -o picrust2_output -p 4 --verbose 
add_descriptions.py -i picrust2_output/pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o picrust2_output/pathways_out/path_abun_unstrat.tsv.descriptions.gz




