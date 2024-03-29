#!/bin/bash
#$ -S /bin/bash
#$ -N diversity_analysis 
#$ -q bioinfo.q
#$ -V 
#$ -cwd 
#$ -notify 
#$ -pe shared 40


#set -e 

source activate qiime2-2020.6
export PERL5LIB='/gpfs0/bioinfo/users/obayomi/miniconda3/envs/qiime2-2020.6/lib/site_perl/5.26.2/x86_64-linux-thread-multi'
export TEMPDIR='/gpfs0/bioinfo/users/obayomi/hinuman_analysis/18S_illumina/tmp/' TMPDIR='/gpfs0/bioinfo/users/obayomi/hinuman_analysis/18S_illumina/tmp/'


TREE=('06.make_tree/dada2' '06.make_tree/dada2' '06.make_tree/dada2' '06.make_tree/deblur' '06.make_tree/deblur' '06.make_tree/dada2' '06.make_tree/dada2' '06.make_tree/dada2' '06.make_tree/deblur' '06.make_tree/deblur')

DEPTH=(437 405 606 122 736 437 405 606 1531 736)

FEATURE_TABLE_DIR=('05.filter_table/dada2' '05.filter_table/dada2'  '05.filter_table/dada2' '05.filter_table/deblur/' '05.filter_table/deblur/' '05.filter_table/dada2/outdoors' '05.filter_table/dada2/outdoors' '05.filter_table/dada2/outdoors' '05.filter_table/deblur/outdoors' '05.filter_table/deblur/outdoors')

PREFIX=('pe' 'se' 'pear-joined' 'se' 'pear-joined' 'pe'  'se' 'pear-joined' 'se' 'pear-joined')

METADATA=('00.mapping/combined.tsv' '00.mapping/combined.tsv' '00.mapping/combined.tsv' '00.mapping/combined.tsv' '00.mapping/combined.tsv' '00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv')

OUT_DIR=('08.core_diversity/dada2' '08.core_diversity/dada2' '08.core_diversity/dada2' '08.core_diversity/deblur' '08.core_diversity/deblur' '08.core_diversity/dada2/outdoors' '08.core_diversity/dada2/outdoors' '08.core_diversity/dada2/outdoors'  '08.core_diversity/deblur/outdoors' '08.core_diversity/deblur/outdoors')

METADATA_COLUMN=('treatment' 'treatment' 'treatment' 'treatment' 'treatment' 'treatment' 'treatment' 'treatment' 'treatment' 'treatment')

# Perform core diversity analysis
parallel --jobs 0 --link qiime diversity core-metrics-phylogenetic \
				--p-sampling-depth {1} \
				--i-table {2}/{3}-filtered_table.qza \
				--i-phylogeny {6}/{3}-rooted-tree.qza \
				--m-metadata-file {4} \
				--output-dir {5}/{3}-diversity-{1} \
				::: ${DEPTH[*]} ::: ${FEATURE_TABLE_DIR[*]} ::: ${PREFIX[*]} ::: ${METADATA[*]} ::: ${OUT_DIR[*]} ::: ${TREE[*]}

# Alpha rarefaction curves show taxon accumulation as a function of sequence depth
parallel --jobs 0 --link qiime diversity alpha-rarefaction \
				--i-table {2}/{3}-filtered_table.qza \
				--p-max-depth {1} \
				--o-visualization {5}/{3}-diversity-{1}/alpha_rarefaction.qzv \
				--m-metadata-file {4} \
				--i-phylogeny {6}/{3}-rooted-tree.qza \
				::: ${DEPTH[*]} ::: ${FEATURE_TABLE_DIR[*]} ::: ${PREFIX[*]} ::: ${METADATA[*]} ::: ${OUT_DIR[*]} ::: ${TREE[*]}

# Alpha Diversity - statistics
# Test for between-group differences

# Faith's phylogenetic diversity
parallel --jobs 0 --link qiime diversity alpha-group-significance \
				--i-alpha-diversity {3}/{1}-diversity-{4}/faith_pd_vector.qza \
				--m-metadata-file {2} \
				--o-visualization {3}/{1}-diversity-{4}/alpha_faith_pd_significance.qzv \
				::: ${PREFIX[*]} ::: ${METADATA[*]} ::: ${OUT_DIR[*]} ::: ${DEPTH[*]}
			
# Shannon Diversity
parallel --jobs 0 --link qiime diversity alpha-group-significance \
				--i-alpha-diversity {3}/{1}-diversity-{4}/shannon_vector.qza \
				--m-metadata-file {2} \
				--o-visualization {3}/{1}-diversity-{4}/alpha_shannon_significance.qzv \
				::: ${PREFIX[*]} ::: ${METADATA[*]} ::: ${OUT_DIR[*]} ::: ${DEPTH[*]}

# Eveness or Chao1
parallel --jobs 0 --link qiime diversity alpha-group-significance \
                                --i-alpha-diversity {3}/{1}-diversity-{4}/evenness_vector.qza \
                                --m-metadata-file {2} \
                                --o-visualization {3}/{1}-diversity-{4}/alpha_evenness_significance.qzv \
                                ::: ${PREFIX[*]} ::: ${METADATA[*]} ::: ${OUT_DIR[*]} ::: ${DEPTH[*]}

# Observed features / OTUs / ASVs
parallel --jobs 0 --link qiime diversity alpha-group-significance \
                                --i-alpha-diversity {3}/{1}-diversity-{4}/observed_features_vector.qza \
                                --m-metadata-file {2} \
                                --o-visualization {3}/{1}-diversity-{4}/alpha_observed_features_significance.qzv \
                                ::: ${PREFIX[*]} ::: ${METADATA[*]} ::: ${OUT_DIR[*]} ::: ${DEPTH[*]}







# Beta Diversity - statistics
# Bray Curtis
parallel --jobs 0 --link qiime diversity beta-group-significance \
				--i-distance-matrix {3}/{1}-diversity-{5}/bray_curtis_distance_matrix.qza \
				--m-metadata-file {2} \
				--m-metadata-column {4} \
				--p-pairwise \
				--o-visualization {3}/{1}-diversity-{5}/beta_bray_curtis_{4}_significance.qzv \
				::: ${PREFIX[*]} ::: ${METADATA[*]} ::: ${OUT_DIR[*]} ::: ${METADATA_COLUMN[*]} ::: ${DEPTH[*]}

# Jaccard
parallel --jobs 0 --link qiime diversity beta-group-significance \
                                --i-distance-matrix {3}/{1}-diversity-{5}/jaccard_distance_matrix.qza \
                                --m-metadata-file {2} \
                                --m-metadata-column {4} \
                                --p-pairwise \
                                --o-visualization {3}/{1}-diversity-{5}/beta_jaccard_{4}_significance.qzv \
                                ::: ${PREFIX[*]} ::: ${METADATA[*]} ::: ${OUT_DIR[*]} ::: ${METADATA_COLUMN[*]} ::: ${DEPTH[*]}

# Unweighted Unifrac
parallel --jobs 0 --link qiime diversity beta-group-significance \
                                --i-distance-matrix {3}/{1}-diversity-{5}/unweighted_unifrac_distance_matrix.qza \
                                --m-metadata-file {2} \
                                --m-metadata-column {4} \
                                --p-pairwise \
                                --o-visualization {3}/{1}-diversity-{5}/beta_unweighted_unifrac_{4}_significance.qzv \
                                ::: ${PREFIX[*]} ::: ${METADATA[*]} ::: ${OUT_DIR[*]} ::: ${METADATA_COLUMN[*]} ::: ${DEPTH[*]}
# weighted Unifrac
parallel --jobs 0 --link qiime diversity beta-group-significance \
                                --i-distance-matrix {3}/{1}-diversity-{5}/weighted_unifrac_distance_matrix.qza \
                                --m-metadata-file {2} \
                                --m-metadata-column {4} \
                                --p-pairwise \
                                --o-visualization {3}/{1}-diversity-{5}/beta_weighted_unifrac_{4}_significance.qzv \
                                ::: ${PREFIX[*]} ::: ${METADATA[*]} ::: ${OUT_DIR[*]} ::: ${METADATA_COLUMN[*]} ::: ${DEPTH[*]}

