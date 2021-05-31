#! /bin/bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N mock_otus
#$ -pe shared 40


set -e
# modify the environment for Qiime
export PATH="/storage16/app/bioinfo/blast-2.2.22/bin:/fastspace/bioinfo_apps/cdbfasta:/fastspace/bioinfo_apps/microbiomeutil-r20110519/ChimeraSlayer:/fastspace/bioinfo_apps/qiime/usr/local/bin:/fastspace/bioinfo_apps/qiime/local/bin:/bin:/usr/bin" PYTHONPATH="/fastspace/bioinfo_apps/qiime/usr/local" HDF5_DIR="/fastspace/bioinfo_apps/qiime/local/hdf5" RDP_JAR_PATH="/fastspace/bioinfo_apps/qiime/local/rdp_classifier_2.2/rdp_classifier-2.2.jar"QIIME_CONFIG_FP="/fastspace/bioinfo_apps/qiime/qiime_config_biores"

#actual command necessay for rdp to function properly
HDF5_DIR="/fastspace/bioinfo_apps/qiime/local/hdf5" RDP_JAR_PATH="/fastspace/bioinfo_apps/qiime/local/rdp_classifier_2.2/rdp_classifier-2.2.jar"

# Pick open reference otus using sortmena for closed ref otu picking and sumaclust for denovo otu picking 
pick_open_reference_otus.py \
		-i /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/04.chimera_check/non_chimeric_seqs.fna \
		-r /gpfs0/bioinfo/users/obayomi/databases/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta \
		-o /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus/ \
		-p /gpfs0/bioinfo/users/obayomi/parameter_files/combined_parameter.txt \
		-m sortmerna_sumaclust \
		-f \
		--suppress_taxonomy_assignment

# Align seqs
#align_seqs.py -i /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus/rep_set.fna \
#	       -o /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//pynast_aligned_seqs \
#		--template_fp /gpfs0/bioinfoi/users/obayomi/databases/SILVA_128_QIIME_release/core_alignment/core_alignment_SILVA128.fna \
#		--alignment_method pynast

# Filter alignment command 
#filter_alignment.py \
#		-o /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//pynast_aligned_seqs \
#		-i /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//pynast_aligned_seqs/rep_set_aligned.fasta

# Build phylogenetic tree command 
#make_phylogeny.py \
#		-i /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta \
#		-o /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//rep_set.tre 


#assign taxonomy using rdp classifier
assign_taxonomy.py \
		-o /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//rdp_assigned_taxonomy \
		-i /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//rep_set.fna \
		--confidence 0.5 \
		--id_to_taxonomy_fp /gpfs0/bioinfo/users/obayomi/databases//SILVA_128_QIIME_release/taxonomy/16S_only/97/consensus_taxonomy_7_levels.txt \
		--assignment_method rdp \
		-r /gpfs0/bioinfo/users/obayomi/databases//SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta \
		--rdp_max_memory 50000

#make otu table
make_otu_table.py \
		-i /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//final_otu_map_mc2.txt  \
		-t /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//rdp_assigned_taxonomy/rep_set_tax_assignments.txt  \
		-o /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//otu_table_MC2.biom \
		-m /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/00.mapping_file_validation/corrected/combined_corrected.txt

# Remove non-target OTUs
filter_taxa_from_otu_table.py \
		-i /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//otu_table_MC2.biom \
		-o /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//otu_table_MC2_taxa_filtered.biom \
		-n D_4__Mitochondria,D_2__Chloroplast,D_0__Archaea,D_0__Eukaryota,Unassigned

# Filter otutable for specific analysis by samples in a mapping files
# water bacteiria (NOT APPLICABLE)
#filter_samples_from_otu_table.py \
#		-i /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//otu_table_MC2_taxa_filtered.biom  \
#		-o /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//combined_otu_table.biom \
#		--sample_id_fp  /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/00.mapping_file_validation/corrected/combined_corrected.txt


# Remove rare otus
filter_otus_from_otu_table.py \
		-i /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//otu_table_MC2_taxa_filtered.biom \
		-o /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//otu_table_MC2_filtered.biom \
		--min_count_fraction 0.00005

# Summarize biom table
biom summarize-table \
		-i /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//otu_table_MC2_filtered.biom \
		-o /gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_16S_illumina/05.open_ref_otus//otu_table_MC2_filtered_biom.summary
 
