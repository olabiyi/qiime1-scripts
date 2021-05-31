#!/bin/bash
#$ -S /bin/bash
#$ -N Filter_samples 
#$ -q bioinfo.q
#$ -V 
#$ -cwd 
#$ -notify 
#$ -pe shared 10

set -e 

source activate qiime2-2020.6
export PERL5LIB='/gpfs0/bioinfo/users/obayomi/miniconda3/envs/qiime2-2020.6/lib/site_perl/5.26.2/x86_64-linux-thread-multi'
export TEMPDIR='/gpfs0/bioinfo/users/obayomi/hinuman_analysis/18S_illumina/tmp/' TMPDIR='/gpfs0/bioinfo/users/obayomi/hinuman_analysis/18S_illumina/tmp/'


OUT_PREFIX=('05.filter_table/dada2/outdoors/pe' '05.filter_table/dada2/outdoors/se' '05.filter_table/dada2/outdoors/pear-joined' '05.filter_table/deblur/outdoors/se' '05.filter_table/deblur/outdoors/pear-joined')

METADATA=('00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv' '00.mapping/outdoors.tsv')

COMBINED_TABLE=('05.filter_table/dada2/pe' '05.filter_table/dada2/se' '05.filter_table/dada2/pear-joined' '05.filter_table/deblur/se' '05.filter_table/deblur/pear-joined')

parallel --jobs 0 --link qiime feature-table filter-samples \
				--i-table {1}-taxa_filtered_table.qza \
				--m-metadata-file {2} \
				--o-filtered-table {3}-taxa_filtered_table.qza ::: ${COMBINED_TABLE[*]} ::: ${METADATA[*]} ::: ${OUT_PREFIX[*]}


parallel --jobs 0 --link qiime feature-table summarize \
                		--i-table {}-taxa_filtered_table.qza \
                		--o-visualization {}-taxa_filtered_table.qzv ::: ${OUT_PREFIX[*]}
