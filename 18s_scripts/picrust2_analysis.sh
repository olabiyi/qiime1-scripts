#!/bin/bash
#$ -S /bin/bash
#$ -N function_analysis
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -notify
#$ -pe shared 40


set -e

# Edit the headers of rep_set.fna to contain only OTU names
#sed -i -E 's/(>.+) .+$/\1/g' rep_set.fna

# make annotation directory
#mkdir /gpfs0/bioinfo/users/obayomi/hinuman_analysis/18S_illumina/11.function_annotation/

##################### Export and rename feature tables and representative sequences from qiime2 artifact

###### copy and rename the artifacts to the function annotation directory
# Feature Tables
cp /gpfs0/bioinfo/users/obayomi/hinuman_analysis/18S_illumina/05.filter_table/dada2/pear-joined-taxa_filtered_table.qza \
  11.function_annotation/
 
# Representative sequences
cp /gpfs0/bioinfo/users/obayomi/hinuman_analysis/18S_illumina/03.dada_denoise/pear-joined-representative_sequences.qza \
 11.function_annotation/ 


cd  /gpfs0/bioinfo/users/obayomi/hinuman_analysis/18S_illumina/11.function_annotation/
source activate qiime2-2020.6
qiime tools export --input-path pear-joined-taxa_filtered_table.qza --output-path ./
mv feature-table.biom pear-joined-feature-table.biom


qiime tools export --input-path pear-joined-representative_sequences.qza --output-path ./
mv dna-sequences.fasta pear-joined-rep_set.fna



source activate picrust2
PREFIX=("pear-joined")
REP_SET=("rep_set.fna")
FEATURE_TABLE=("feature-table.biom")


function run_picrust(){

	local PREFIX=$1
	local REP_SET=$2
	local FEATURE_TABLE=$3
	# Run PICRUST2 pipeline 
	picrust2_pipeline.py \
		-s ${PREFIX}-${REP_SET} \
		-i ${PREFIX}-${FEATURE_TABLE} \
		-o ${PREFIX}-picrust2_out_pipeline \
		-p 40

	# Annotate you enzymes / pathways by adding a description column
	add_descriptions.py -i ${PREFIX}-picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                	    -o ${PREFIX}-picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

	add_descriptions.py -i ${PREFIX}-picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
        	            -o ${PREFIX}-picrust2_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv.gz

	add_descriptions.py -i ${PREFIX}-picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
        	            -o ${PREFIX}-picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

	# Unzip the prediction files
	gunzip ${PREFIX}-picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
	gunzip ${PREFIX}-picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
	gunzip ${PREFIX}-picrust2_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv.gz

	gunzip ${PREFIX}-picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz
	gunzip ${PREFIX}-picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz
	gunzip ${PREFIX}-picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv.gz

	# Convert to biom
	biom convert \
		-i ${PREFIX}-picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv \
		-o ${PREFIX}-picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.biom \
		--table-type="OTU table" \
		--to-hdf5

	biom convert \
		-i ${PREFIX}-picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv \
		-o ${PREFIX}-picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.biom \
		--table-type="OTU table" \
		--to-hdf5

	biom convert \
		-i ${PREFIX}-picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv \
		-o ${PREFIX}-picrust2_out_pipeline/pathways_out/path_abun_unstrat.biom \
		--table-type="OTU table" \
		--to-hdf5

}


export -f run_picrust

parallel --jobs 0 --link run_picrust {1} {2} {3} ::: ${PREFIX[*]} ::: ${REP_SET[*]} ::: ${FEATURE_TABLE[*]}
