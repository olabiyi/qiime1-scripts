#! /bin/bash
#$ -S /bin/bash
#$ -q bioinfo.q
#$ -V
#$ -cwd
#$ -N pathogens_search
#$ -pe shared 24



###### This script runs steps 2-3 in the 'how_to_run_pathogenAnalysis.txt' fiile in one go on qsub ######


# Source your .bashrc file - it should look like the lines below
#this is just to make available the necessary scripts for the analysis
#pathogen analysis scripts
#export PATH=/path/to/the/folder/containing/the/scripts/for/pathogen_analysis:$PATH
#qiime
#export PATH=/path/to/your/qiime/usr/local/bin/:$PATH
#NCBI blast
#export PATH=/path/to/your/ncbi-blast-2.10.1+/bin/:$PATH

source /gpfs0/bioinfo/users/obayomi/.bashrc
mapping_file=/gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_18S_illumina/00.mapping_file_validation/corrected/combined_corrected.txt
# Get the list of samples in a QIIME mapping file
#by getting the required field, removing the first line then
#sorting it, getting the unique sample names and saving them in a file called sample_ids.txt
cut -f1 ${mapping_file} |sed '1d'| sort -u  > sample_ids.txt 

# Assign the samples to a variable
samples=$(cat sample_ids.txt)

# Convert to samples variable to an array
samples=($samples)


# colName assumes there is a column in your mapping file called sample, therefore create it
# Enter the following necessary  variables MANUALLY
search=pathogens
pathogens=/gpfs0/bioinfo/users/obayomi/hinuman_analysis/16s_pathogen_analysis/human_protozoa_list.txt
database=/gpfs0/bioinfo/users/obayomi/databases/REF_SEQ_18s/18S_all
fasta=/gpfs0/bioinfo/users/obayomi/hinuman_analysis/q1_18S_illumina/07.pathogen_analysis/protist/pathogens.fasta
colName=sample
results_together=results_together
ID=9

# Leave these set variables as they are
seqs=_$search.fasta
blast_out=_organisms.blast
list=_$search.list
organism_list=_organisms.list
unique=_unique
identity=_identity
final=final_
analyses=_Analysis

#loop for as many samples
#run a blast search on every treatment / sample
for treatment in ${samples[*]}; do
   # Remove the directory if it already exists and make a new directory for the treatment and then change to the directory
   rm -rf $treatment$analyses
   mkdir $treatment$analyses
   cd $treatment$analyses
    
    # Extract  sequences that belong to the treatment
   extract_seqs_by_sample_id.py -i $fasta -o $treatment$seqs \
   -m ${mapping_file} -s "$colName:$treatment"

   echo "blasting... $treatment"
   # Blast the treatments's sequences against any blast formatted database, retrieving only the best hit
   blastn -db $database -query  $treatment$seqs -outfmt 0 \
   -num_descriptions 1 -num_alignments 1 -out $treatment$blast_out

   # Get the list of organisms from the blast output with their identities
    grep ">"  $treatment$blast_out > $treatment$organism_list

   # Get the list of the organisms' identity scores
   grep -E "Identities" $treatment$blast_out > $treatment$identity
   # Paste the organism list and Identities together
   paste -d" " $treatment$organism_list $treatment$identity > $treatment$organism_list$identity 
   # Get the list of organisms with identity scores greater than ID or 100%
        if [ $ID -eq 10 ]; then
             grep -E "(100%)" $treatment$organism_list$identity > $treatment$identity$organism_list
        else 
             grep -E "(9[$ID-9]%)|(100%)" $treatment$organism_list$identity > $treatment$identity$organism_list
        fi
   # Create a file that will contain the list of what you are looking for
   touch $treatment$identity$list
   # Search for what you are looking for within the organisms list
   grep -wFf $pathogens  $treatment$identity$organism_list > $treatment$identity$list
   ## Get the pathogen list and parse it for readability i.e. without the sequence assertion number and removing anything 16S or 18S
   cut -d" " -f2-6  $treatment$identity$list |sed -e 's/1[68]S//g'|sort > $final$treatment$list

   # Get the unique list
   sort -u $final$treatment$list > $treatment$unique$list
   # Return to the parent directory
   cd ..
done

#####   Organize the output files ##################
# Remove the folders containing the combined together if they already exist
rm -rf $results_together
rm -rf unique

# Create a directory that will contain the final results 
mkdir $results_together
mkdir unique
# Put all the final and unique results files in one folder 
find . -type f -name "final_*" -exec mv -vuni '{}' "$results_together" ";" 2> error.log 
find . -type f -name "*unique*" -exec mv -vuni '{}' "unique/" ";" 2>> error.log

# Change the extensions of the final and unique lists from .list to .txt so that they may be viewed 
# with a text editor outside unix such as notepad
rename .list .txt $results_together/*
rename .list .txt unique/*
