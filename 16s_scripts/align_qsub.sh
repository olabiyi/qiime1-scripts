#! /bin/bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N wash_align
#$ -pe shared 40

# modify the environment for Qiime
export PATH="/storage16/app/bioinfo/blast-2.2.22/bin:/fastspace/bioinfo_apps/cdbfasta:/fastspace/bioinfo_apps/microbiomeutil-r20110519/ChimeraSlayer:/fastspace/bioinfo_apps/qiime/usr/local/bin:/fastspace/bioinfo_apps/qiime/local/bin:/bin:/usr/bin" PYTHONPATH="/fastspace/bioinfo_apps/qiime/usr/local" HDF5_DIR="/fastspace/bioinfo_apps/qiime/local/hdf5" RDP_JAR_PATH="/fastspace/bioinfo_apps/qiime/local/rdp_classifier_2.2/rdp_classifier-2.2.jar"QIIME_CONFIG_FP="/fastspace/bioinfo_apps/qiime/qiime_config_biores"

#align sequences
align_seqs.py -i  /gpfs0/biores/users/gilloro/Biyi/community_analysis/fruit_wash/05.open_ref_otus/rep_set.fna -t /gpfs0/biores/users/gilloro/Biyi/SILVA_128_QIIME_release/core_alignment/core_alignment_SILVA128.fna -o /gpfs0/biores/users/gilloro/Biyi/community_analysis/fruit_wash/05.open_ref_otus/pynast_aligned_seqs/
