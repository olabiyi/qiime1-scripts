#! /bin/bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N usearch_chimera_check
#$ -pe shared 40

set -e

# modify the environment for Qiime
export PATH="/storage16/app/bioinfo/blast-2.2.22/bin:/fastspace/bioinfo_apps/cdbfasta:/fastspace/bioinfo_apps/microbiomeutil-r20110519/ChimeraSlayer:/fastspace/bioinfo_apps/qiime/usr/local/bin:/fastspace/bioinfo_apps/qiime/local/bin:/bin:/usr/bin" PYTHONPATH="/fastspace/bioinfo_apps/qiime/usr/local" HDF5_DIR="/fastspace/bioinfo_apps/qiime/local/hdf5" RDP_JAR_PATH="/fastspace/bioinfo_apps/qiime/local/rdp_classifier_2.2/rdp_classifier-2.2.jar"QIIME_CONFIG_FP="/fastspace/bioinfo_apps/qiime/qiime_config_biores" 
#PERL5LIB='/gpfs0/bioinfo/users/obayomi/miniconda3/envs/qiime2-2020.6/lib/site_perl/5.26.2/x86_64-linux-thread-multi'


export DATABASE="/gpfs0/bioinfo/users/obayomi/databases/SILVA_128_QIIME_release/rep_set/rep_set_18S_only/97/97_otus_18S.fasta"
export DIR="/gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_18S_illumina/04.chimera_check/"
export SPLIT_FASTA="${DIR}/split_samples/"
export SEQS_FNA="/gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_18S_illumina/03.demultiplexed/seqs.fna"


#chimera checking with usearch61

# split the seqs.fna file into per sample fna files - 
# this step is necessary because usearch fails when you have a lot of samples
split_sequence_file_on_sample_ids.py \
	-i ${SEQS_FNA} \
	-o ${SPLIT_FASTA}


FILES=($(find "${SPLIT_FASTA}" -type f -name "*.fasta"))


function chimera_check(){

	local file=$1
        # Remove anything before and including split_samples/
	local sample=${file#*split_samples/}
	sample=${sample/.fasta/}
	
	mkdir ${SPLIT_FASTA}/${sample}

	identify_chimeric_seqs.py \
	-m usearch61 \
	-i ${file} \
	-r ${DATABASE} \
	-o ${SPLIT_FASTA}/${sample}

	cat ${SPLIT_FASTA}/${sample}/chimeras.txt \
	>> ${DIR}/CHIMERA.txt

	cat ${SPLIT_FASTA}/${sample}/non_chimeras.txt \
	>> ${DIR}/NON_CHIMERA.txt

}

export -f chimera_check

#source activate qiime2-2020.6
#PERL5LIB='/gpfs0/bioinfo/users/obayomi/miniconda3/envs/qiime2-2020.6/lib/site_perl/5.26.2/x86_64-linux-thread-multi'

/gpfs0/bioinfo/users/obayomi/bin/parallel --jobs 0 chimera_check {} ::: ${FILES[*]}

#remove unnessary files
rm -rf  ${SPLIT_FASTA}

#remove chimeric sequences fom seqs.fna
filter_fasta.py \
	-s ${DIR}/CHIMERA.txt \
	-f ${SEQS_FNA} \
	-o ${DIR}/non_chimeric_seqs.fna \
	-n
