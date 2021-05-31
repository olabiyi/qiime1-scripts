#! /bin/bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N wash_usearch_chimera_check
#$ -pe shared 40

# modify the environment for Qiime
export PATH="/storage16/app/bioinfo/blast-2.2.22/bin:/fastspace/bioinfo_apps/cdbfasta:/fastspace/bioinfo_apps/microbiomeutil-r20110519/ChimeraSlayer:/fastspace/bioinfo_apps/qiime/usr/local/bin:/fastspace/bioinfo_apps/qiime/local/bin:/bin:/usr/bin" PYTHONPATH="/fastspace/bioinfo_apps/qiime/usr/local" HDF5_DIR="/fastspace/bioinfo_apps/qiime/local/hdf5" RDP_JAR_PATH="/fastspace/bioinfo_apps/qiime/local/rdp_classifier_2.2/rdp_classifier-2.2.jar"QIIME_CONFIG_FP="/fastspace/bioinfo_apps/qiime/qiime_config_biores"

#identify chimera
#identify_chimeric_seqs.py -m ChimeraSlayer -i /gpfs0/biores/users/gilloro/Biyi/community_analysis/water/05.open_ref_otus/bac_open_ref_otus/pynast_aligned_seqs/rep_set_aligned.fasta -a /gpfs0/biores/users/gilloro/Biyi/SILVA_128_QIIME_release/core_alignment/core_alignment_SILVA128.fna -o /gpfs0/biores/users/gilloro/Biyi/community_analysis/water/05.open_ref_otus/bac_open_ref_otus/chimera_checking/chimeric_seqs.txt

#identify chimera using vsearch

#source /storage/users/gilloro/.bashrc

#run chimera check
#vsearch --uchime_denovo /gpfs0/biores/users/gilloro/Biyi/community_analysis/water/03.demultiplexed/bacteria_demultiplexed/seqs.derep.fna -strand plus --chimeras /gpfs0/biores/users/gilloro/Biyi/community_analysis/water/05.open_ref_otus/bac_open_ref_otus/chimera_checking/chimeric_seqs.fasta --nonchimeras /gpfs0/biores/users/gilloro/Biyi/community_analysis/water/05.open_ref_otus/bac_open_ref_otus/chimera_checking/non_chimeric_seqs.fasta --log /gpfs0/biores/users/gilloro/Biyi/community_analysis/water/05.open_ref_otus/bac_open_ref_otus/chimera_checking/denovo.log




#loop through every sample while ignoring out files from qsub
for file in $(ls -I "wash*"); do 
identify_chimeric_seqs.py -m usearch61 -i $file -r /gpfs0/biores/users/gilloro/Biyi/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta -o usearch61_chimera_checking/
cat usearch61_chimera_checking/chimeras.txt >> usearch61_chimera_checking/CHIMERA.txt
cat usearch61_chimera_checking/non_chimeras.txt >> usearch61_chimera_checking/NON_CHIMERA.txt
done

#remove unnessary files
rm -rf usearch61_chimera_checking/*.fasta; rm -rf usearch61_chimera_checking/*.uc; rm -rf usearch61_chimera_checking/*.log; rm -rf usearch61_chimera_checking/*.uchime;
rm -rf usearch61_chimera_checking/chimeras.txt ; rm -rf usearch61_chimera_checking/non_chimeras.txt
