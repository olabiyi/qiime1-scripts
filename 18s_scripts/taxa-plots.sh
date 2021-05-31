#!/bin/bash
#$ -S /bin/bash
#$ -N make_bar_plots 
#$ -q bioinfo.q
#$ -V 
#$ -cwd 
#$ -notify 
#$ -pe shared 40

set -e 

source activate qiime2-2020.6
export PERL5LIB='/gpfs0/bioinfo/users/obayomi/miniconda3/envs/qiime2-2020.6/lib/site_perl/5.26.2/x86_64-linux-thread-multi'

export TEMPDIR='/gpfs0/bioinfo/users/obayomi/hinuman_analysis/18S_illumina/tmp/' TMPDIR='/gpfs0/bioinfo/users/obayomi/hinuman_analysis/18S_illumina/tmp/'


TAXONOMY_DIR=('04.assign_taxonomy/dada2' '04.assign_taxonomy/dada2' '04.assign_taxonomy/dada2' '04.assign_taxonomy/deblur' '04.assign_taxonomy/deblur' '04.assign_taxonomy/dada2' '04.assign_taxonomy/dada2' '04.assign_taxonomy/dada2' '04.assign_taxonomy/deblur' '04.assign_taxonomy/deblur')

FEATURE_TABLE_DIR=('05.filter_table/dada2' '05.filter_table/dada2' '05.filter_table/dada2' '05.filter_table/deblur' '05.filter_table/deblur' '05.filter_table/dada2/outdoors' '05.filter_table/dada2/outdoors' '05.filter_table/dada2/outdoors' '05.filter_table/deblur/outdoors' '05.filter_table/deblur/outdoors')

PREFIX=('pe' 'se' 'pear-joined' 'se' 'pear-joined' 'pe' 'se' 'pear-joined' 'se' 'pear-joined')

METADATA=('00.mapping/combined.tsv' '00.mapping/combined.tsv' '00.mapping/combined.tsv' '00.mapping/combined.tsv' '00.mapping/combined.tsv' '00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv')

METADATA_COLUMN=('treatment' 'treatment' 'treatment' 'treatment' 'treatment' 'treatment' 'treatment' 'treatment' 'treatment' 'treatment')

MODE=('combined' 'combined' 'combined' 'combined' 'combined' 'outdoors' 'outdoors' 'outdoors' 'outdoors' 'outdoors')

GROUP_METADATA=('00.mapping/combined-treatment.tsv' '00.mapping/combined-treatment.tsv' '00.mapping/combined-treatment.tsv' '00.mapping/combined-treatment.tsv' '00.mapping/combined-treatment.tsv' '00.mapping/outdoors-treatment.tsv' '00.mapping/outdoors-treatment.tsv' '00.mapping/outdoors-treatment.tsv' '00.mapping/outdoors-treatment.tsv' '00.mapping/outdoors-treatment.tsv')

PLOT_DIR=('07.make_taxa_plots/dada2' '07.make_taxa_plots/dada2' '07.make_taxa_plots/dada2' '07.make_taxa_plots/deblur' '07.make_taxa_plots/deblur' '07.make_taxa_plots/dada2/outdoors' '07.make_taxa_plots/dada2/outdoors' '07.make_taxa_plots/dada2/outdoors' '07.make_taxa_plots/deblur/outdoors' '07.make_taxa_plots/deblur/outdoors')




# Sample bar plots
parallel --jobs 0 --link qiime taxa barplot \
				--i-table {1}/{2}-filtered_table.qza \
				--i-taxonomy {6}/{2}-taxonomy.qza \
				--m-metadata-file {3} \
				--o-visualization  {4}/{5}-samples-{2}-bar-plots.qzv \
				::: ${FEATURE_TABLE_DIR[*]} ::: ${PREFIX[*]} ::: ${METADATA[*]} ::: ${PLOT_DIR[*]} ::: ${MODE[*]} ::: ${TAXONOMY_DIR[*]}

# Taxa bar plots of metadata group - here by treatment
# group feature table (*-filtered_table.qza)  by metadata column of interest
parallel --jobs 0 --link qiime feature-table group \
				--i-table  {1}/{2}-filtered_table.qza \
				--p-axis sample \
				--m-metadata-file {3} \
				--m-metadata-column {4} \
				--p-mode sum \
				--o-grouped-table {1}/{5}-{4}-{2}-filtered_table.qza \
				::: ${FEATURE_TABLE_DIR[*]} ::: ${PREFIX[*]} ::: ${METADATA[*]} ::: ${METADATA_COLUMN[*]} ::: ${MODE[*]}

# Make bar plot of group table here by treatment
# Make sure you create a new metadata with the group level names as sample-id
parallel --jobs 0 --link qiime taxa barplot \
				--i-table  {1}/{5}-{4}-{2}-filtered_table.qza \
				--i-taxonomy {6}/{2}-taxonomy.qza \
				--m-metadata-file {3} \
				--o-visualization {7}/{5}-{4}-{2}-bar-plots.qzv \
				::: ${FEATURE_TABLE_DIR[*]} ::: ${PREFIX[*]} ::: ${GROUP_METADATA[*]} ::: ${METADATA_COLUMN[*]} ::: ${MODE[*]} ::: ${TAXONOMY_DIR[*]} ::: ${PLOT_DIR[*]}

