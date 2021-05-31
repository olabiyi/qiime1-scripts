#!/bin/bash
#$ -S /bin/bash
#$ -N Filter_features 
#$ -q bioinfo.q
#$ -V 
#$ -cwd 
#$ -notify 
#$ -pe shared 10

set -e 

source activate qiime2-2020.6
export PERL5LIB='/gpfs0/bioinfo/users/obayomi/miniconda3/envs/qiime2-2020.6/lib/site_perl/5.26.2/x86_64-linux-thread-multi'
export TEMPDIR='/gpfs0/bioinfo/users/obayomi/hinuman_analysis/18S_illumina/tmp/' TMPDIR='/gpfs0/bioinfo/users/obayomi/hinuman_analysis/18S_illumina/tmp/'
export EXCLUDE="Bacteria,Fungi,Chytridiomycota,Basidiomycota,Metazoa,Rotifera,Gastrotricha,Nematozoa,Embryophyta,Spermatophyta,Asterales,Brassicales,Caryophyllales,Cupressales,Fabales,Malpighiales,Pinales,Rosales,Solanales,Arecales,Asparagales,Poales,Capsicum,Jatropha,Bryophyta,Tracheophyta"

IN_PREFIX=('03.dada_denoise/pe' '03.dada_denoise/se' '03.dada_denoise/pear-joined' '03.deblur_denoise/se' '03.deblur_denoise/pear-joined')

# For Combined table i.e the original table with indoors, outdoors and mock tables combined
#OUT_PREFIX=('05.filter_table/dada2/pe' '05.filter_table/dada2/se' '05.filter_table/dada2/pear-joined' '05.filter_table/deblur/se' '05.filter_table/deblur/pear-joined')
TAXONOMY_PREFIX=('04.assign_taxonomy/dada2/pe' '04.assign_taxonomy/dada2/se' '04.assign_taxonomy/dada2/pear-joined' '04.assign_taxonomy/deblur/se' '04.assign_taxonomy/deblur/pear-joined')
#TOTAL_SEQUENCES=(83640 173756 210131 15865 9230) multiply each number by 0.00005 to get the minimum number for filtering rare otus below
#MIN_FREQUENCY=(4 9 11 1 1)

# For the tables that have been split by metadata
OUT_PREFIX=('05.filter_table/dada2/outdoors/pe' '05.filter_table/dada2/outdoors/se' '05.filter_table/dada2/outdoors/pear-joined' '05.filter_table/deblur/outdoors/se' '05.filter_table/deblur/outdoors/pear-joined')
#TOTAL_SEQUENCES=(67107 149549 183754 12315 6379)
MIN_FREQUENCY=(3 7 9 1 1)

REMOVE_RARE_FEATURES="true" #"false"

function filter_table(){
	
	local in_prefix=$1
	local out_prefix=$2
	local taxonomy_prefix=$3

	# Remove singletons
	qiime feature-table filter-features \
			--i-table ${in_prefix}-table.qza \
			--p-min-frequency 2 \
			--o-filtered-table ${out_prefix}-noSingleton_filtered_table.qza

	qiime feature-table summarize \
			--i-table ${out_prefix}-noSingleton_filtered_table.qza \
			--o-visualization ${out_prefix}-noSingleton_filtered_table.qzv


	# Remove unassigned, archaea, eukaryota, chloroplast and mitochondria taxa
	qiime taxa filter-table \
		--i-table ${out_prefix}-noSingleton_filtered_table.qza \
		--i-taxonomy ${taxonomy_prefix}-taxonomy.qza \
		--p-exclude ${EXCLUDE} \
		--o-filtered-table ${out_prefix}-taxa_filtered_table.qza

	# To figure out the total number of sequences ("Total freqency") here equals ${TOTAL_SEQUENCES} e.g. 8,053,326
	qiime feature-table summarize \
		--i-table ${out_prefix}-taxa_filtered_table.qza \
		--o-visualization ${out_prefix}-taxa_filtered_table.qzv

}

if [ "${REMOVE_RARE_FEATURES}" == "false" ]; then

	export -f filter_table
	parallel --jobs 0 --link filter_table {1} {2} {3} ::: ${IN_PREFIX[*]} ::: ${OUT_PREFIX[*]} ::: ${TAXONOMY_PREFIX[*]}

else
	##### Removing rare otus / features with abundance less the 0.005%
	parallel --jobs 0 --link qiime feature-table filter-features \
  				--i-table {1}-taxa_filtered_table.qza \
  				--p-min-frequency {2} \
  				--o-filtered-table {1}-filtered_table.qza ::: ${OUT_PREFIX[*]} ::: ${MIN_FREQUENCY[*]}

	parallel --jobs 0 --link qiime feature-table summarize \
                --i-table {}-filtered_table.qza \
                --o-visualization {}-filtered_table.qzv ::: ${OUT_PREFIX[*]}

fi
