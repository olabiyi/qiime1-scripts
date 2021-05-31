#! /bin/bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N demultiplex
#$ -pe shared 10

set -e

# modify the environment for Qiime
export PATH="/storage16/app/bioinfo/blast-2.2.22/bin:/fastspace/bioinfo_apps/cdbfasta:/fastspace/bioinfo_apps/microbiomeutil-r20110519/ChimeraSlayer:/fastspace/bioinfo_apps/qiime/usr/local/bin:/fastspace/bioinfo_apps/qiime/local/bin:/bin:/usr/bin" PYTHONPATH="/fastspace/bioinfo_apps/qiime/usr/local" HDF5_DIR="/fastspace/bioinfo_apps/qiime/local/hdf5" RDP_JAR_PATH="/fastspace/bioinfo_apps/qiime/local/rdp_classifier_2.2/rdp_classifier-2.2.jar"QIIME_CONFIG_FP="/fastspace/bioinfo_apps/qiime/qiime_config_biores"

#actual command
HDF5_DIR="/fastspace/bioinfo_apps/qiime/local/hdf5" RDP_JAR_PATH="/fastspace/bioinfo_apps/qiime/local/rdp_classifier_2.2/rdp_classifier-2.2.jar"



# concatenate the joined sequences in one seqs.fna file - You can add the -w flag for debugging
multiple_split_libraries_fastq.py \
	-i /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/02.join/ \
	-o /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/03.demultiplexed/ \
	-p /gpfs0/bioinfo/users/obayomi/parameter_files/combined_parameter.txt

#get the statistics for split libray log (OPTIONAL)
#this step is optional but highly recommended
Rscript /gpfs0/bioinfo/users/obayomi/split_library_log_parser/parse_split_library_log.R 
	--input /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/03.demultiplexed/split_library_log.txt 
	--output /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/03.demultiplexed/split_library_log.tsv
