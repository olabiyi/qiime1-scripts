#! /bin/bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N coreDiversity
#$ -pe shared 40

# modify the environment for Qiime
export PATH="/storage16/app/bioinfo/blast-2.2.22/bin:/fastspace/bioinfo_apps/cdbfasta:/fastspace/bioinfo_apps/microbiomeutil-r20110519/ChimeraSlayer:/fastspace/bioinfo_apps/qiime/usr/local/bin:/fastspace/bioinfo_apps/qiime/local/bin:/bin:/usr/bin" PYTHONPATH="/fastspace/bioinfo_apps/qiime/usr/local" HDF5_DIR="/fastspace/bioinfo_apps/qiime/local/hdf5" RDP_JAR_PATH="/fastspace/bioinfo_apps/qiime/local/rdp_classifier_2.2/rdp_classifier-2.2.jar"QIIME_CONFIG_FP="/fastspace/bioinfo_apps/qiime/qiime_config_biores"


DEPTH=(2770 15744)
# Perform corediversity analysis - alpha and beta diversity and relative abundance bar-plots
core_diversity_analyses.py \
		-i /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus/otu_table_MC2_filtered.biom \
		-o /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/06.core_diversity/diversity-15744/ \
		--mapping_fp /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/00.mapping_file_validation/corrected/mock.tsv_corrected.txt \
		--sampling_depth 15744 \
		--parameter_fp /gpfs0/bioinfo/users/obayomi/parameter_files/mock_parameter.txt \
		--tree_fp /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus/rep_set.tre \
		--categories treatment
