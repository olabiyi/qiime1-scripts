#!/bin/bash
#$ -S /bin/bash
#$ -N diff_abundance
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -notify
#$ -pe shared 40


set -e

source activate qiime2-2020.6
export PERL5LIB='/gpfs0/bioinfo/users/obayomi/miniconda3/envs/qiime2-2020.6/lib/site_perl/5.26.2/x86_64-linux-thread-multi'
export TEMPDIR='/gpfs0/bioinfo/users/obayomi/hinuman_analysis/18S_illumina/tmp/' TMPDIR='/gpfs0/bioinfo/users/obayomi/hinuman_analysis/18S_illumina/tmp/'

TAXON_LEVELS=(2 3 4 5 6)

FEATURE_TABLE_DIR=('05.filter_table/dada2' '05.filter_table/dada2' '05.filter_table/dada2'  '05.filter_table/deblur/' '05.filter_table/deblur/' '05.filter_table/dada2/outdoors' '05.filter_table/dada2/outdoors' '05.filter_table/dada2/outdoors'  '05.filter_table/deblur/outdoors' '05.filter_table/deblur/outdoors')

PREFIX=('pe' 'se' 'pear-joined' 'se' 'pear-joined' 'pe' 'se' 'pear-joined' 'se' 'pear-joined')

METADATA=('00.mapping/combined.tsv' '00.mapping/combined.tsv' '00.mapping/combined.tsv' '00.mapping/combined.tsv' '00.mapping/combined.tsv' '00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv')

OUT_DIR=('09.differential_abundance/dada2' '09.differential_abundance/dada2' '09.differential_abundance/dada2' '09.differential_abundance/deblur' '09.differential_abundance/deblur' '09.differential_abundance/dada2/outdoors' '09.differential_abundance/dada2/outdoors' '09.differential_abundance/dada2/outdoors' '09.differential_abundance/deblur/outdoors' '09.differential_abundance/deblur/outdoors')

METADATA_COLUMN=('treatment' 'treatment' 'treatment' 'treatment' 'treatment' 'treatment' 'treatment' 'treatment' 'treatment' 'treatment')

TAXONOMY_DIR=('04.assign_taxonomy/dada2' '04.assign_taxonomy/dada2' '04.assign_taxonomy/dada2'  '04.assign_taxonomy/deblur' '04.assign_taxonomy/deblur' '04.assign_taxonomy/dada2' '04.assign_taxonomy/dada2' '04.assign_taxonomy/dada2'  '04.assign_taxonomy/deblur' '04.assign_taxonomy/deblur')



# Differntial abundance testing using ANCOM
# At ASV level
# Add pseudocount to ASV table because ANCOM can't deal with zero counts
parallel --jobs 0 --link qiime composition add-pseudocount \
				--i-table {1}/{2}-filtered_table.qza \
				--o-composition-table {3}/{2}-composition-table.qza \
				::: ${FEATURE_TABLE_DIR[*]} ::: ${PREFIX[*]} ::: ${OUT_DIR[*]}  
  
# Apply ANCOM to identify ASV/OTUs that differ in abundance
parallel --jobs 0 --link qiime composition ancom  \
				--i-table {3}/{1}-composition-table.qza \
				--m-metadata-file {2} \
				--m-metadata-column {4} \
				--o-visualization {3}/{1}-{4}-ancom.qzv \
				::: ${PREFIX[*]} ::: ${METADATA[*]} ::: ${OUT_DIR[*]} ::: ${METADATA_COLUMN[*]}



for TAXON_LEVEL in ${TAXON_LEVELS[*]}; do

	# At specific taxonomy level - here at the genus level level 6 (L6)
	# 1.  Collapse feauture table at a taxonomy level of interest
	parallel --jobs 0 --link qiime taxa collapse \
				--i-table {1}/{2}-filtered_table.qza \
				--i-taxonomy {4}/{2}-taxonomy.qza \
				--p-level ${TAXON_LEVEL} \
				--o-collapsed-table {3}/{2}-L${TAXON_LEVEL}-filtered_table.qza \
				::: ${FEATURE_TABLE_DIR[*]} ::: ${PREFIX[*]} ::: ${OUT_DIR[*]} ::: ${TAXONOMY_DIR[*]}
				

	# 2. Add pseudocount to ASV table because ANCOM can't deal with zero counts
	parallel --jobs 0 --link qiime composition add-pseudocount \
				--i-table {2}/{1}-L${TAXON_LEVEL}-filtered_table.qza \
				--o-composition-table {2}/{1}-L${TAXON_LEVEL}-composition-table.qza \
				::: ${PREFIX[*]} ::: ${OUT_DIR[*]}
  
	# 3. Apply ANCOM to identify ASV/OTUs that differ in abundance
	parallel --jobs 0 --link qiime composition ancom  \
				--i-table {3}/{1}-L${TAXON_LEVEL}-composition-table.qza \
				--m-metadata-file {2} \
				--m-metadata-column {4} \
				--o-visualization {3}/{1}-L${TAXON_LEVEL}-{4}-ancom.qzv \
			       ::: ${PREFIX[*]} ::: ${METADATA[*]} ::: ${OUT_DIR[*]} ::: ${METADATA_COLUMN[*]}

done

