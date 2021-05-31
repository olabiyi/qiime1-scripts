#!/bin/bash
#$ -S /bin/bash
#$ -N denoize_dada2 
#$ -q bioinfo.q
#$ -V 
#$ -cwd 
#$ -notify 
#$ -pe shared 40

set -e

source activate qiime2-2020.6
PAIRED="false"
#PAIRED="true"
TRIM_LEFT=0
TRIM_RIGHT=0
TRUNC_LENGTH=400
#TRUNC_LENGTH=260
TRUNC_LENGTH_LEFT=260
TRUNC_LENGTH_RIGHT=240
maxE_f=4
maxE_r=7

#PREFIX="se"
#PREFIX="pe"
PREFIX="pear-joined"

if [ "${PAIRED}" != "true" ]; then
	echo "running dada single"
	# Denoise, truncate and assign ASVs
	qiime dada2 denoise-single \
  		--i-demultiplexed-seqs 01.import/${PREFIX}-reads.qza \
  		--p-trim-left ${TRIM_LEFT} \
  		--p-trunc-len ${TRUNC_LENGTH} \
		--p-max-ee ${maxE_f} \
  		--o-representative-sequences 03.dada_denoise/${PREFIX}-representative_sequences.qza \
  		--o-table 03.dada_denoise/${PREFIX}-table.qza \
  		--o-denoising-stats  03.dada_denoise/${PREFIX}-denoise_stats.qza


	qiime feature-table summarize \
		--i-table 03.dada_denoise/${PREFIX}-table.qza \
		--o-visualization 03.dada_denoise/${PREFIX}-table_summary.qzv


	qiime feature-table tabulate-seqs \
  		--i-data 03.dada_denoise/${PREFIX}-representative_sequences.qza \
  		--o-visualization 03.dada_denoise/${PREFIX}-representative_sequences.qzv



	qiime metadata tabulate \
		--m-input-file 03.dada_denoise/${PREFIX}-denoise_stats.qza \
 		--o-visualization 03.dada_denoise/${PREFIX}-denoise_stats.qzv

else
		echo "running dada paired"
	qiime dada2 denoise-paired \
		--i-demultiplexed-seqs 01.import/${PREFIX}-reads.qza \
		--o-table 03.dada_denoise/${PREFIX}-table.qza \
		--o-representative-sequences 03.dada_denoise/${PREFIX}-representative_sequences.qza \
		--o-denoising-stats 03.dada_denoise/${PREFIX}-denoise_stats.qza \
		--p-trunc-len-f ${TRUNC_LENGTH_LEFT} \
		--p-trunc-len-r ${TRUNC_LENGTH_RIGHT} \
		--p-trim-left-f ${TRIM_LEFT} \
		--p-trim-left-r ${TRIM_RIGHT} \
		--p-max-ee-f ${maxE_f} \
		--p-max-ee-r ${maxE_r} \
		--p-n-threads 30


	# This visualization shows us the sequences per sample spread - to determine minimum number for rarefaction
	# and sequences per feature (OTU or ASV)
	qiime feature-table summarize \
		--i-table 03.dada_denoise/${PREFIX}-table.qza \
		--o-visualization 03.dada_denoise/${PREFIX}-table_summary.qzv


	qiime feature-table tabulate-seqs \
  		--i-data 03.dada_denoise/${PREFIX}-representative_sequences.qza \
  		--o-visualization 03.dada_denoise/${PREFIX}-representative_sequences.qzv


	qiime metadata tabulate \
		--m-input-file 03.dada_denoise/${PREFIX}-denoise_stats.qza \
		--o-visualization 03.dada_denoise/${PREFIX}-denoise_stats.qzv

fi
