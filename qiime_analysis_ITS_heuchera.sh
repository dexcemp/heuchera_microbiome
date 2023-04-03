
# activate conda virtual environment
conda activate qiime2-2020.8 

# Make a manifest file of this format
######################
sample-id	forward-absolute-filepath	reverse-absolute-filepath
A40	$PWD/test_seqs/A40_P1.fastq.gz	$PWD/test_seqs/A40_P2.fastq.gz
######################


## Build UNITE database for fungal ITS (this was already done on BotBot) -- these steps take a LONG time -- about 2 hours
## If there are no classifications below kingdom for a custom database, consider subsampling to one sequence for the host as below.
#conda activate qiime2-2020.8 
#cd UNITE
#awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' developer/sh_refs_qiime_ver8_99_02.02.2019_dev.fasta  | tr -d ' ' > developer/sh_refs_qiime_ver8_99_02.02.2019_dev_uppercase.fasta
#
#cat developer/sh_refs_qiime_ver8_99_02.02.2019_dev_uppercase.fasta Heuchera_ITS_database.one.fasta > unite_heuchera_99.fasta
#qiime tools import --type FeatureData[Sequence] --input-path unite_heuchera_99.fasta --output-path unite_heuchera_99.qza
#cat developer/sh_taxonomy_qiime_ver8_99_02.02.2019_dev.txt Heuchera_ITS_database.one.csv > unite_heuchera_99.txt
#qiime tools import --type FeatureData[Taxonomy] --input-path unite_heuchera_99.txt --output-path unite_heuchera_99-tax.qza --input-format HeaderlessTSVTaxonomyFormat
## Skip the reference trimming, unlike 16S, per Qiime documentation
#qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads unite_heuchera_99.qza --i-reference-taxonomy unite_heuchera_99-tax.qza --o-classifier unite_heuchera_99-classifier.qza

# import sequences
# Qiime must have gzipped input files in a folder
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.txt --output-path import.qza --input-format PairedEndFastqManifestPhred33V2

# Trimming is to remove primers, set to ITS primers ITS1FI2/ITS2
qiime dada2 denoise-paired --i-demultiplexed-seqs import.qza --p-trim-left-f 18 --p-trim-left-r 20 --p-trunc-len-f 250 --p-trunc-len-r 200 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza

# Taxonomic analysis
qiime feature-classifier classify-sklearn --i-classifier ./UNITE/unite_heuchera_99-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

# Remove chloroplast and mitochondrial 16S
qiime taxa filter-table --p-exclude k__Viridiplantae --i-table table.qza --i-taxonomy taxonomy.qza --o-filtered-table table.hostremoved.qza

qiime taxa filter-seqs --p-exclude k__Viridiplantae --i-sequences rep-seqs.qza --i-taxonomy taxonomy.qza --o-filtered-sequences rep-seqs.hostremoved.qza



# Taxa barplots, no host
qiime taxa barplot --i-table table.hostremoved.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.tsv --o-visualization taxa-bar-plots.qzv




# Run standard feature tables
qiime feature-table summarize --i-table table.hostremoved.qza --o-visualization table.qzv --m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs --i-data rep-seqs.hostremoved.qza --o-visualization rep-seqs.qzv

qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv

qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.hostremoved.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

# Diversity statistics and significance
# p-sampling-depth: change per experiment; look across your samples for feature counts -- choose a number as high as possible but lower thant the lowest sample
# Will not work for too few samples
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.hostremoved.qza --p-sampling-depth 11 --output-dir core-metrics-results --m-metadata-file sample-metadata.tsv 
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/faith_pd_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization core-metrics-results/faith-pd-group-significance.qzv
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/evenness_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization core-metrics-results/evenness-group-significance.qzv
# Due to the longer time needed for beta diversity significance, a specific metadata column is specified to reduce computational time
qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column host_species --o-visualization core-metrics-results/unweighted-unifrac-host_species-significance.qzv --p-pairwise
qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column host_subsection --o-visualization core-metrics-results/unweighted-unifrac-host_subsection-significance.qzv --p-pairwise

# Ordination plots
qiime emperor plot --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza --m-metadata-file sample-metadata.tsv --p-custom-axes host_species --o-visualization core-metrics-results/unweighted-unifrac-emperor-host_species.qzv
qiime emperor plot --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza --m-metadata-file sample-metadata.tsv --p-custom-axes host_species --o-visualization core-metrics-results/bray-curtis-emperor-host_species.qzv

# Rarefaction plots
# p-sampling-depth: change per experiment; look across your samples for feature counts -- choose a number as high as possible but lower thant the lowest sample
qiime diversity alpha-rarefaction --i-table table.hostremoved.qza --i-phylogeny rooted-tree.qza --p-max-depth 11 --m-metadata-file sample-metadata.tsv --o-visualization alpha-rarefaction.qzv

# ANCOM -- differential abundance analysis
# Helps find which genera, families, etc. differ across metadata categories
qiime feature-table filter-samples --i-table table.hostremoved.qza --m-metadata-file sample-metadata.tsv --p-where "[host_species]='Heuchera_americana_americana'" --o-filtered-table host_species-table.qza
# Imputation method for zero-abundance samples
qiime composition add-pseudocount --i-table host_species-table.qza --o-composition-table comp-host_species-table.qza
qiime composition ancom --i-table comp-host_species-table.qza --m-metadata-file sample-metadata.tsv --m-metadata-column host_species --o-visualization ancom-host_species.qzv

# ANCOM with taxonomic level collapsing -- here level 6 genus
qiime taxa collapse --i-table host_species-table.qza --i-taxonomy taxonomy.qza --p-level 6 --o-collapsed-table host_species-table-l6.qza
qiime composition add-pseudocount --i-table host_species-table-l6.qza --o-composition-table comp-host_species-table-l6.qza
qiime composition ancom --i-table comp-host_species-table-l6.qza --m-metadata-file sample-metadata.tsv --m-metadata-column host_species --o-visualization l6-ancom-host_species.qzv

# adonis
qiime diversity adonis --i-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --p-formula "host_species+host_subsection" --p-permutations 10000 --p-n-jobs 4 --o-visualization ./core-metrics-results/adonis.unweighted_unifrac.qzv



# Analyse how much was thrown away as host DNA using new taxonomic barplots
qiime feature-classifier classify-sklearn --i-classifier ./UNITE/unite_heuchera_99-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.withhost.qza
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.withhost.qzv
qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.withhost.qza --m-metadata-file sample-metadata.tsv --o-visualization taxa-bar-plots.withhost.qzv


#######
# Aggregate taxon barplots to host species

# don't forget the .tsv file mentioned below -- it's a simple spreadsheet

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

