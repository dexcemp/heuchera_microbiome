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
qiime greengenes2 filter-features --i-feature-table table.qza --i-reference 2022.10.taxonomy.asv.nwk.qza --o-filtered-feature-table table.gg2.qza
qiime greengenes2 taxonomy-from-table --i-reference-taxonomy 2022.10.taxonomy.asv.nwk.qza --i-table table.gg2.qza --o-classification taxonomy.qza

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

# Remove chloroplast and mitochondrial 16S
qiime taxa filter-table --p-exclude c__Chloroplast,f__mitochondria --i-table table.qza --i-taxonomy taxonomy.qza --o-filtered-table table.organellesremoved.qza
qiime taxa filter-seqs --p-exclude c__Chloroplast,f__mitochondria --i-sequences rep-seqs.qza --i-taxonomy taxonomy.qza --o-filtered-sequences rep-seqs.organellesremoved.qza

# Taxa barplots, no organelles
qiime taxa barplot --i-table table.organellesremoved.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.tsv --o-visualization taxa-bar-plots.qzv

# Rarefaction plots
# p-sampling-depth: change per experiment; look across your samples for feature counts -- choose a number as high as possible but lower thant the lowest sample
qiime diversity alpha-rarefaction --i-table table.organellesremoved.qza --i-phylogeny rooted-tree.qza --p-max-depth 50 --m-metadata-file sample-metadata.tsv --o-visualization alpha-rarefaction.qzv
# with host
qiime diversity alpha-rarefaction --i-table table.qza --i-phylogeny rooted-tree.qza --p-max-depth 150 --m-metadata-file sample-metadata.tsv --o-visualization alpha-rarefaction.withhost150.qzv
qiime diversity alpha-rarefaction --i-table table.qza --i-phylogeny rooted-tree.qza --p-max-depth 500 --m-metadata-file sample-metadata.tsv --o-visualization alpha-rarefaction.withhost500.qzv
qiime diversity alpha-rarefaction --i-table table.qza --i-phylogeny rooted-tree.qza --p-max-depth 1000 --m-metadata-file sample-metadata.tsv --o-visualization alpha-rarefaction.withhost1000.qzv

