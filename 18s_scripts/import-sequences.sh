#!/bin/bash
#$ -S /bin/bash
#$ -N import_reads 
#$ -q bioinfo.q
#$ -V 
#$ -cwd 
#$ -notify 
#$ -pe shared 10

set -e

source activate qiime2-2020.6
export PERL5LIB='/gpfs0/bioinfo/users/obayomi/miniconda3/envs/qiime2-2020.6/lib/site_perl/5.26.2/x86_64-linux-thread-multi'
export IMPORT_DIR='01.import' QC_DIR='02.QC'


PREFIX=("pe" "se" "pear-joined")
READ_TYPE=("paired-end" "single-end" "pear-joined")
MANIFEST_FILES=("sequence_data/pe-MANIFEST" "sequence_data/se-MANIFEST" "stitched_reads/joined-MANIFEST")

function import_reads(){

	local READ_TYPE=$1
	local MANIFEST=$2
	local PREFIX=$3

	if [ "${READ_TYPE}" == "paired-end" ]; then

		# import pe-sequnces to qiime
		qiime tools import \
			--type 'SampleData[PairedEndSequencesWithQuality]' \
			--input-path ${MANIFEST} \
			--output-path ${IMPORT_DIR}/${PREFIX}-reads.qza \
			--input-format PairedEndFastqManifestPhred33


	else

		# import se-sequnces or joined reads to qiime
		qiime tools import \
			--type 'SampleData[SequencesWithQuality]' \
			--input-path  ${MANIFEST} \
			--output-path ${IMPORT_DIR}/${PREFIX}-reads.qza \
			--input-format SingleEndFastqManifestPhred33

	fi

	# Demultiplex and View reads quality
	# Analyze quality scores of 10000 random samples
	# paired end
	qiime demux summarize \
		--p-n 10000 \
		--i-data ${IMPORT_DIR}/${PREFIX}-reads.qza \
		--o-visualization ${QC_DIR}/${PREFIX}-qual_viz.qzv


}


export -f import_reads

parallel --jobs 0 --link import_reads {} {} {} ::: ${READ_TYPE[*]} ::: ${MANIFEST_FILES[*]} ::: ${PREFIX[*]}
