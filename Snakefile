from os import path, getcwd
import pandas as pd


# ------- A pipeline to perform microbiome and function profiling using QIIME 1.9.1 ----------- #
# Run on your local computer like so:
# snakemake -npr --cores 10 --keep-going --rerun-incomplete --restart-times 3
configfile: "config/config.yaml"
container: config['containers']['qiime']
onsuccess:
    print("Workflow completed without any error")
onerror:
    print("An error occurred")

rule all:
    input:
        "02.Count_seqs_pre_trim/seqs_stat.txt",
        "03.QC/multiqc_report.html",
        "07.SummarizeQC_post_adaptor_trim/multiqc_report.html",
        "08.Count_seqs_post_trim/seqs_stat.txt",
        "11.Merge_reads/read_counts.tsv",
        "12.Demultiplexed/seqs_count.txt",
        "19.Summarize_table/otu_table_rare_taxa_filtered.biom.summary",
        "20.Diversity_analysis_{depth}/index.html".format(depth=config['parameters']['core_diversity']['depth']),
        "21.Function_annotation/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv",
        "23.Filter_otu2seq_map/pathogens_map.txt",
        "24.Filter_fasta/pathogens.fasta" 
        


# --------- Quality check, trimming and adapter removal  --------------- #

rule Count_seqs_pre_trim:
    input: expand(["01.raw_data/{sample}/{sample}_R1.fastq.gz", \
                   "01.raw_data/{sample}/{sample}_R2.fastq.gz"], \
                    sample=config['samples'])
    output: "02.Count_seqs_pre_trim/seqs_stat.txt"
    container: config['containers']['seqkit']
    shell:
        """
        # Get the stats on the sequences using seqkit
        seqkit stats {input} > temp.txt

         # Sort the sequence statistics
         (sed -n '1p' temp.txt; awk 'NR>1{{print}}' temp.txt | \
           sort -V -k1,1) > {output} \
           && rm temp.txt
        """


rule QC_pre_trim:
    input:
        forward="01.raw_data/{sample}/{sample}_R1.fastq.gz",
        rev="01.raw_data/{sample}/{sample}_R2.fastq.gz"
    output:
        forward_html="03.QC/{sample}/{sample}_R1_fastqc.html",
        rev_html="03.QC/{sample}/{sample}_R2_fastqc.html"
    container: config['containers']['fastqc']
    params:
        out_dir=lambda w, output: path.dirname(output[0]),
        threads=5
    log: "logs/QC_pre_trim/{sample}/{sample}.log"
    threads: 5
    shell:
        "fastqc --outdir {params.out_dir}/ "
        "--threads {params.threads} {input.forward} {input.rev}"

rule SummarizeQC_pre_trim:
    input:
        forward_html=expand("03.QC/{sample}/{sample}_R1_fastqc.html",
                            sample=config['samples']),
        rev_html=expand("03.QC/{sample}/{sample}_R2_fastqc.html",
                            sample=config['samples'])
    output: "03.QC/multiqc_report.html"
    log: "logs/SummarizeQC_pre_trim/multiqc.log"
    container: config['containers']['multiqc']
    params:
        out_dir=lambda w, output: path.dirname(output[0])
    threads: 1
    shell:
        "multiqc --interactive -f {params.out_dir} -o {params.out_dir}"


adaptors=config['parameters']['trimmomatic']['adaptors']
min_length=config['parameters']['trimmomatic']['min_len']
rule Trimmomatic_trim_adaptors:
    input:
        forward="01.raw_data/{sample}/{sample}_R1.fastq.gz",
        rev="01.raw_data/{sample}/{sample}_R2.fastq.gz",
    output:
        r1="04.Trimmomatic_trim_adaptors/{sample}/{sample}_R1.fastq.gz",
        r2="04.Trimmomatic_trim_adaptors/{sample}/{sample}_R2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="04.Trimmomatic_trim_adaptors/{sample}/{sample}_R1.unpaired.fastq.gz",
        r2_unpaired="04.Trimmomatic_trim_adaptors/{sample}/{sample}_R2.unpaired.fastq.gz"
    log: "logs/Trimmomatic_trim_adaptors/{sample}/{sample}.log"
    conda: config['containers']['trimmomatic']
    params:
        trimmer="ILLUMINACLIP:{adaptors}:2:30:10"
                " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20"
                " MINLEN:{min_length}".format(adaptors=adaptors,
                                          min_length=min_length)
    threads: 1
    resources:
        mem_mb=1024
    shell:
        "trimmomatic PE "
        "-threads {threads} "
        "{input.forward} {input.rev} "
        "{output.r1} {output.r1_unpaired} "
        "{output.r2} {output.r2_unpaired} "
        "{params.trimmer} > {log} 2>&1 "

# Trim Primers using cutadapt
rule Trim_primers:
    input:
        forward_reads=rules.Trimmomatic_trim_adaptors.output.r1,
        rev_reads=rules.Trimmomatic_trim_adaptors.output.r2
    output:
        forward_reads="05.Trim_primers/{sample}/{sample}_R1.fastq.gz",
        rev_reads="05.Trim_primers/{sample}/{sample}_R2.fastq.gz"
    log: "logs/Trim_primers/{sample}/{sample}.log"
    container: config['containers']['cutadapt']
    params:
        forward_primer=config['parameters']['cutadapt']['forward_primer'],
        rev_primer=config['parameters']['cutadapt']['reverse_primer'],
        minimum_length=config['parameters']['cutadapt']['minimum_length'],
        quality_cutoff=config['parameters']['cutadapt']['quality_cutoff']
    threads: 1
    shell:
        """
        # -a TGGAATTCTCGGGTGCCAAGG sequence here is the small RNA 3' adaptor
        # https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/TruSeq/RNA/Small-RNA/TruSeqSmallRNA.htm
        # -A CTGTCTCTTATACAC sequece here is the nextera transposa sequence 
        cutadapt \
              -g '{params.forward_primer}' \
              -G '{params.rev_primer}' \
              -o {output.forward_reads} \
              -p {output.rev_reads} \
              --minimum-length  {params.minimum_length} \
              --quality-cutoff  {params.quality_cutoff} \
             {input.forward_reads} {input.rev_reads} > {log} 2>&1

        """

# Trim adaptors and quality check using cutadapt and factqc with Trim galore
rule Trim_galore_trim_adaptors:
    input: 
        forward=rules.Trim_primers.output.forward_reads,
        rev=rules.Trim_primers.output.rev_reads
    output:
        forward_reads="06.Trim_galore_trim_adaptors/{sample}/{sample}_R1.fastq.gz",
        rev_reads="06.Trim_galore_trim_adaptors/{sample}/{sample}_R2.fastq.gz",
        forward_html="06.Trim_galore_trim_adaptors/{sample}/{sample}_R1_fastqc.html",
        rev_html="06.Trim_galore_trim_adaptors/{sample}/{sample}_R2_fastqc.html"
    log: "logs/Trim_galore_trim_adaptors/{sample}/{sample}.log"
    threads: 1
    container: config['containers']['trim_galore']
    params:
        out_dir=lambda w, output: path.dirname(output[0])
    shell:
        """ 
         trim_galore \
           -o {params.out_dir} \
           --fastqc  \
           --paired {input.forward} {input.rev} > {log} 2>&1

         #Rename the files
         #Fastq files
        mv {params.out_dir}/{wildcards.sample}_R1_val_1.fq.gz {params.out_dir}/{wildcards.sample}_R1.fastq.gz
        mv {params.out_dir}/{wildcards.sample}_R2_val_2.fq.gz {params.out_dir}/{wildcards.sample}_R2.fastq.gz
         #HTML files
        mv {params.out_dir}/{wildcards.sample}_R1_val_1_fastqc.html {params.out_dir}/{wildcards.sample}_R1_fastqc.html
        mv {params.out_dir}/{wildcards.sample}_R2_val_2_fastqc.html {params.out_dir}/{wildcards.sample}_R2_fastqc.html
         #Zip files
        mv {params.out_dir}/{wildcards.sample}_R1_val_1_fastqc.zip {params.out_dir}/{wildcards.sample}_R1_fastqc.zip
        mv {params.out_dir}/{wildcards.sample}_R2_val_2_fastqc.zip {params.out_dir}/{wildcards.sample}_R2_fastqc.zip
        """


# Aggregate and summarize the quality check using multiqc
rule SummarizeQC_post_adaptor_trim:
    input:
        forward_html=expand("06.Trim_galore_trim_adaptors/{sample}/{sample}_R1_fastqc.html",
                            sample=config['samples']),
        rev_html=expand("06.Trim_galore_trim_adaptors/{sample}/{sample}_R2_fastqc.html",
                            sample=config['samples'])
    output: "07.SummarizeQC_post_adaptor_trim/multiqc_report.html"
    log: "logs/SummarizeQC_post_adaptor_trim/multiqc.log"
    container: config['containers']['multiqc']
    params:
        in_dir=lambda w, input: path.dirname(input[0]).split("/")[0],
        out_dir=lambda w, output: path.dirname(output[0])
    threads: 1
    shell:
        """
          multiqc \
              --interactive \
              -f {params.in_dir} \
              -o {params.out_dir}  > {log} 2>&1
        """


rule Count_seqs_post_trim:
    input: expand(["06.Trim_galore_trim_adaptors/{sample}/{sample}_R1.fastq.gz", \
                   "06.Trim_galore_trim_adaptors/{sample}/{sample}_R2.fastq.gz"], \
                   sample=config['samples'])
    output: "08.Count_seqs_post_trim/seqs_stat.txt"
    container: config['containers']['seqkit']
    shell:
        """
        # Get the stats on the sequences using seqkit
        seqkit stats {input} > temp.txt

         # Sort the sequence statistics
         (sed -n '1p' temp.txt; awk 'NR>1{{print}}' temp.txt | \
           sort -V -k1,1) > {output} \
           && rm temp.txt
        """
#### ---------------------------- Qiime Begins ---------------------------------------------------###

# Rename the mapping file so that it will be easy to work with
rule Rename_mapping_file:
    input: config['mapping_file']
    output: "09.Rename_mapping_file/mapping_file"
    log: "logs/Rename_mapping_file/Rename_mapping_file.log"
    threads: 1
    shell:
        "cat {input} > {output}"

rule Validate_mapping_file:
    input: rules.Rename_mapping_file.output 
    output: "10.Validate_mapping_file/corrected/mapping_file_corrected.txt"
    log: "logs/Validate_mapping_file/Validate_mapping_file.log"
    threads: 1
    params:
        out_dir=lambda w, output: path.dirname(output[0]),
        basename=lambda w, output: path.basename(output[0])
    shell:
        """
        validate_mapping_file.py \
           -m {input} \
           -o {params.out_dir}
        """

rule Merge_reads:
    input:
        forward="06.Trim_galore_trim_adaptors/{sample}/{sample}_R1.fastq.gz",
        rev="06.Trim_galore_trim_adaptors/{sample}/{sample}_R2.fastq.gz"
    output: "11.Merge_reads/{sample}/{sample}.fastq"
    log: "logs/Merge_reads/{sample}/{sample}.log"
    threads: 5
    container: config['containers']['pear']
    params: 
        out_dir=lambda w, output: path.dirname(output[0]),
        min=config['parameters']['pear']['min_assembly'],
        max=config['parameters']['pear']['max_assembly'],
        min_trim=config['parameters']['pear']['min_trim'],
        threads=config['parameters']['pear']['threads']
    shell:
        """    
        [ -d {params.out_dir} ] ||  mkdir -p {params.out_dir}
        # Merge reads then delete unnecessary files
        pear \
            -f {input.forward} \
            -r {input.rev} \
            -j {params.threads} \
            -o {params.out_dir}/{wildcards.sample} \
            -m {params.max} \
            -n {params.min} \
            -t {params.min_trim} > {log} 2>&1
        
       rm -rf \
            {params.out_dir}/{wildcards.sample}.discarded.fastq \
            {params.out_dir}/{wildcards.sample}.unassembled.forward.fastq \
            {params.out_dir}/{wildcards.sample}.unassembled.reverse.fastq 
            mv {params.out_dir}/{wildcards.sample}.assembled.fastq {params.out_dir}/{wildcards.sample}.fastq
         
       # gzip to save memory
         
       #gzip {params.out_dir}/{wildcards.sample}.fastq
        """

rule CountSeqs:
    input: expand("11.Merge_reads/{sample}/{sample}.fastq", \
                  sample=config['samples'])
    output: "11.Merge_reads/read_counts.tsv"
    container: config['containers']['seqkit']
    log: "logs/CountSeqs/CountSeqs.log"
    threads: 2
    shell:
        """
        seqkit stats {input} > {output} 2> {log}
        """

# Concatenate and demultiplex reads into one seqs.fna file
rule Demultiplex_reads:
    input:
        reads=expand(rules.Merge_reads.output, sample=config['samples']),
        parameter_file=ancient(config['parameter_file']),
        mapping_file=rules.Validate_mapping_file.output
    output: 
        seqs="12.Demultiplexed/seqs.fna",
        log="12.Demultiplexed/split_library_log.txt"
    log: "logs/Demultiplex_reads/Demultiplex_reads.log"
    threads: 5
    #container: config['containers']['qiime']
    params:
        out_dir=lambda w, output: path.dirname(output[0]),
        reads_dir=lambda w, input: path.dirname(input.reads[0]).split("/")[0]
    shell:
        """
        multiple_split_libraries_fastq.py --sampleid_indicator '.fastq' \
             -i {params.reads_dir} -o {params.out_dir} \
             -p {input.parameter_file}
        """
        
# Parse the split library log file in order to retrieve statistics 
# about the demultiplexing process
#rule Split_library_statistics:
#    input: rules.Demultiplex_reads.output.log
#    output: "06.demultiplexed/split_library_log.tsv"
#    log: "logs/Split_library_statistics/Split_library_statistics.log"
#    threads: 1
#    params:
#        conda_activate=config['QIIME_ENV'],
#        program=config['programs_path']['split_library_log_parser']
#    shell:
#        """
#        Rscript {params.program} --input {input} --output {output}
#        """

rule Count_sequences:
    input: rules.Demultiplex_reads.output.seqs
    output: "12.Demultiplexed/seqs_count.txt"
    log: "logs/Count_sequences/Count_sequences.log"
    threads: 1
    shell:
        """
        # Count all the sequences after split_labrary's quality filtering and concatenation
        count_seqs.py -i {input} > {output}
        """

# ----------------- Chimera detection and filtering using usearch61 ------------------------ #

# Because Usearch61 is is not free and open source
# In order to use the free version we need to split
# our seqs.fna file in smaller chunks
# then indentify chimeras in each chunk

rule Split_seqs_file:
    input: rules.Demultiplex_reads.output.seqs
    output: expand("13.Split_sequences/split_samples/{sample}.fasta", \
                   sample=config['samples'])
    log: "logs/Split_seqs_file/Split_seqs_file.log"
    threads: 5
    params: 
       out_dir=lambda w, output: path.dirname(output[0])
    shell:
        """
        # Split the seqs.fna file into per sample fna files 
        # this step is necessary because usearch fails when you have a lot of samples
        split_sequence_file_on_sample_ids.py  -i {input} -o  {params.out_dir}
        """

# Denovo and reference based chimera detection using usearch61
rule Identify_chimeras:
    input: 
        fasta="13.Split_sequences/split_samples/{sample}.fasta",
        reference_database=config['database']
    output: "14.Identify_chimeras/chimera_checking/{sample}_chimera.txt"
    log: "logs/Identify_chimeras/{sample}_Identify_chimeras.log"
    threads: 2
    params:
       out_dir=lambda w, input: path.dirname(input.fasta),
       chimera_dir=lambda w, output: path.dirname(output[0]),
       usearch61_dir="/usr/local/bin/"
    shell:
        """
          
        [ -d {params.out_dir}/{wildcards.sample}/ ] || mkdir -p {params.out_dir}/{wildcards.sample}/
        PROJECT_DIR=`pwd`
        cd {params.out_dir}/{wildcards.sample}/ 

          identify_chimeric_seqs.py \
             -m usearch61 \
             -i ${{PROJECT_DIR}}/{input.fasta} \
             -r {input.reference_database} \
             -o ${{PROJECT_DIR}}/{params.out_dir}/{wildcards.sample}/

        [ -d ${{PROJECT_DIR}}/{params.chimera_dir}/ ] || mkdir -p ${{PROJECT_DIR}}/{params.chimera_dir}/
        # Rename the chimera.txt file so that it will be sample specific
        mv ${{PROJECT_DIR}}/{params.out_dir}/{wildcards.sample}/chimeras.txt \
          ${{PROJECT_DIR}}/{params.chimera_dir}/{wildcards.sample}_chimera.txt && \
        rm -rf  ${{PROJECT_DIR}}/{params.out_dir}/{wildcards.sample}/
        """

rule Concatenate_chimeras:
    input: expand("14.Identify_chimeras/chimera_checking/{sample}_chimera.txt", \
                  sample=config['samples'])
    output: "14.Identify_chimeras/chimera_checking/all_chimera.txt"
    log: "logs/Contatenate_chimeras/Contatenate_chimeras.log"
    threads: 1
    shell:
        "cat {input} > {output}"


rule Filter_chimeras:
    input: 
        seqs=rules.Demultiplex_reads.output.seqs,
        chimeras=rules.Concatenate_chimeras.output
    output: "15.Filter_chimeras/non_chimeric_seqs.fna"
    log: "logs/Filter_chimeras/Filter_chimeras.log"
    threads: 5
    shell:
        """
        # Filter out the chimeric sequences fom seqs.fna
        filter_fasta.py \
           -s {input.chimeras} \
           -f {input.seqs} \
           -o {output} \
           -n
        """

# Count the total number of reads after chimera removal
rule Count_non_chimeric_sequences:
    input: rules.Filter_chimeras.output
    output: "15.Filter_chimeras/non_chimeric_seqs_count.txt"
    log: "logs/Count_non_chimeric_sequences/Count_non_chimeric_sequences.log"
    threads: 1
    shell:
        """
        set +u
        {params.conda_activate}
        set -u

        # Count all the sequences after chimera filtering
        count_seqs.py -i {input} > {output}
        """

# -------------------------------- OTU picking ------------------------------- #
# Pick open-reference OTUs using sortmerna and sumaclust algorithms, 
# construct a phylogenetic tree and assign taxonomy
rule Pick_otus:
    input:
        seqs=rules.Filter_chimeras.output,
        reference_database=config['database'],
        parameter_file=config['parameter_file']
    output:
        otu_table="16.Pick_otus/otus/otu_table_mc2_w_tax.biom",
        representative_sequences="16.Pick_otus/otus/new_refseqs.fna",
        otu_map="16.Pick_otus/otus/final_otu_map.txt",
        #tree="16.Pick_otus/otus/rep_set.tre"
    log: "logs/Pick_otus/Pick_otus.log"
    threads: 10
    params:
       out_dir=lambda w, output: path.dirname(output.otu_table),
       method=config['otu_picking_method'],
       method_class="open_ref"
    shell:
        """
        #set +u;source activate qiime1;set -u
        if [ {params.method_class} == "denovo" ];then
        pick_de_novo_otus.py \
            -i {input.seqs} \
            -o {params.out_dir} \
            -p {input.parameter_file} \
            -f

       else

        pick_open_reference_otus.py \
            -i {input.seqs} \
            -r {input.reference_database} \
            -o {params.out_dir} \
            -p {input.parameter_file} \
            -m {params.method} \
            --suppress_align_and_tree \
            -f

        fi
        """   


amplicon = config['amplicon'] # "16S", "18S", "ITS"

# Filter out non target taxa 
rule Filter_non_target_taxa:
    input: rules.Pick_otus.output.otu_table
    output: "17.Filter_non_target_taxa/otu_table_taxa_filtered.biom"
    log: "logs/Filter_non_target_taxa/Filter_non_target_taxa.log"
    threads: 2
    params:
        taxa2filter=config['taxa2filter']
    shell:
        """
        filter_taxa_from_otu_table.py \
            -i {input} \
            -o {output} \
            -n {params.taxa2filter}
        """

# Filter out rare otus that seem to skew our dataset and make it difficult to interpret our results
rule Filter_rare_otus:
    input: rules.Filter_non_target_taxa.output
    output: "18.Filter_rare_otus/otu_table_rare_taxa_filtered.biom"
    log: "logs/Filter_rare_otus/Filter_rare_otus.log"
    threads: 2
    params:
        min_freq=config['parameters']['filter_rare']['minimum_frequency']
    shell:
        """
        # Filter-out the really rare otus "The recommended procedure is to discard those OTUs with a
        # number of sequences less than 0.005% of the total number of sequences" (Navas-Molina et al, 2013)
        filter_otus_from_otu_table.py \
            -i {input} \
            -o {output} \
            --min_count_fraction {params.min_freq}
        """

# Summarize the biom table that you may know the number of sequences for rarefaction 
# - it is usually the lowest number of sequences but this decision
# is yours to make. I rarefied at 1000 sequences per sample but you can go higher or lower as you so please
# always confirm if rarefaction was fine using rarefaction curves

rule Summarize_table:
    input: rules.Filter_rare_otus.output
    output: "19.Summarize_table/otu_table_rare_taxa_filtered.biom.summary"
    log: "logs/Summarize_table/Summarize_table.log"
    threads: 1
    shell:
        """
        biom summarize-table \
            -i {input} \
            -o {output}
        """

rule Core_diversity_analyses:
    input: 
        otu_table=rules.Filter_rare_otus.output,
        mapping_file=rules.Validate_mapping_file.output,
        #tree=rules.Pick_otus.output.tree
    output: "20.Diversity_analysis_{depth}/index.html".format(depth=config['parameters']['core_diversity']['depth'])
    log: "logs/Core_diversity_analyses/Core_diversity_analyses.log"
    threads: 10
    container: config['containers']['qiime_stan']
    params:
        out_dir=lambda w, output: path.dirname(output[0]),
        depth=config['parameters']['core_diversity']['depth'],
        category=config['parameters']['core_diversity']['category'],
        parameter_file=config['parameter_file']
    shell:
        """
        #  #--tree_fp
        if [ -d {params.out_dir} ]; then rm -rf {params.out_dir}; fi 

        set +u;source activate qiime1;set -u && \
        core_diversity_analyses.py \
		    -i  {input.otu_table} \
		    -o  {params.out_dir} \
		    --mapping_fp {input.mapping_file} \
		    --sampling_depth {params.depth} \
		    --parameter_fp {params.parameter_file} \
		    --categories {params.category} \
            --nonphylogenetic_diversity
        """

if config['amplicon'] == "16S":

    # ------------------ Function analysis using Picrust2 -----------------#

    rule Function_annotation:
        input:
            feature_table=rules.Filter_rare_otus.output,
            rep_seqs=rules.Pick_otus.output.representative_sequences
        output:
            ec="21.Function_annotation/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
            ko="21.Function_annotation/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
            pathway="21.Function_annotation/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv.gz",
            #contrib="21.Function_annotation/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_contrib.tsv.gz"
        log: "logs/Function_annotation/Function_annotation.log"
        threads: 10
        params:
            out_dir=lambda w, output: output.ec.split('/')[0] + "/" + output.ec.split('/')[1],
            conda_activate=config['containers']['picrust2']
        shell:
            """
            {params.conda_activate}
            # Remove the temporary output directory if it already exists
            [ -d picrust2_out_pipeline/ ] && rm -rf picrust2_out_pipeline/
        
            # ---- Run picrust2 pipeline for function annotation -------- #
            picrust2_pipeline.py \
                -s {input.rep_seqs} \
                -i {input.feature_table} \
                -o picrust2_out_pipeline/ \
                -p {threads} && \
                mv picrust2_out_pipeline/* {params.out_dir}/ && \
                rmdir picrust2_out_pipeline/
            """

    # Add description to PICRUST2 function annotation tables
    rule Add_description:
        input:
            ec=rules.Function_annotation.output.ec,
            ko=rules.Function_annotation.output.ko,
            pathway=rules.Function_annotation.output.pathway,
            #contrib=rules.Function_annotation.output.contrib
        output:
            ec="21.Function_annotation/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv",
            ko="21.Function_annotation/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv",
            pathway="21.Function_annotation/picrust2_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv",
            #ec_contrib="21.Function_annotation/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_contrib.tsv",
            #ko_contrib="21.Function_annotation/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_contrib.tsv",
            #pathway_contrib="21.Function_annotation/picrust2_out_pipeline/pathways_out/path_abun_contrib.tsv"
        log: "logs/Add_description/Add_description.log"
        threads: 10
        #conda: "picrust2"  #config['containers']['picrust2']
        params:
            outdir="21.Function_annotation/picrust2_out_pipeline/",
            conda_activate=config['containers']['picrust2']
        shell:
            """
            {params.conda_activate}
            # ----- Annotate your enzymes, KOs and pathways by adding a description column ------#
            # EC
            add_descriptions.py -i {input.ec} -m EC -o {output.ec}
            # Metacyc Pathway
            add_descriptions.py -i {input.pathway} -m METACYC -o {output.pathway}
            # KO
            add_descriptions.py -i {input.ko} -m KO -o {output.ko} 
 
            # Unizip the metagenome contribution files - these files describe the micribes contribution the function profiles
            #find {params.outdir} -type f -name "*contrib.tsv.gz" -exec gunzip {{}} \;
            """

else:
    # This is a dummy rule that just creates an empty file
    rule Add_description:
        input:
            feature_table=rules.Filter_rare_otus.output
        output:
            ko="21.Function_annotation/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv"
        log: "logs/Add_description/Add_description.log"
        threads: 1
        params:
            outdir="21.Function_annotation/picrust2_out_pipeline/KO_metagenome_out"
        shell:
            """
                # Create an empty file
                mkdir -p {params.outdir} && touch {output.ko}
            """


# ------------- Potenetial pathogen analysis using 16S rRNA sequences ---------------- #

# Generate pathogen table
rule Filter_taxa:
    input:
        biom_table=rules.Pick_otus.output.otu_table,
        list2search=config['Pathogen_list']
    output:
        biom="22.Filter_taxa/pathogens.biom",
        tsv="22.Filter_taxa/pathogens.tsv"
    threads: 1
    log: "logs/Filter_taxa/Filter_taxa.log"
    shell:
        """
        # filter out the OTUs of what you are looking for from the given OTU table e.g. human pathogens
        # And convert biom formattted OTU table to tsv 
        filter_taxa_from_otu_table.py \
            -i {input.biom_table} \
            -o {output.biom} -p {input.list2search} && \
        biom convert     \
               -i  {output.biom}   \
               -o  {output.tsv}     \
               --to-tsv     \
               --header-key taxonomy     \
               --output-metadata-id "Consensus Lineage" > {log} 2>&1
        """

# Subset OTU to sequence map to contain
# only potential pathogenic otus
rule Filter_otu2seq_map:
    input:
        tsv_table=rules.Filter_taxa.output.tsv,
        otu_map=rules.Pick_otus.output.otu_map
    output: "23.Filter_otu2seq_map/pathogens_map.txt"
    threads: 10
    log: "logs/Filter_otu2seq_map/Filter_otu2seq_map.log"
    conda: config['containers']['bioinfo']
    shell:
        """
        # Get the pathogenic OTU names
        OTUS=($(cut -f1 {input.tsv_table} | sed -e 1,2d))
        parallel -j {threads} "grep -wE '{{}}' {input.otu_map}  >> {output}" ::: ${{OTUS[*]}} \
            > {log} 2>&1
        """

# Get the pathogens fasta sequences
rule Filter_fasta:
    input:
        seqs=rules.Filter_chimeras.output,
        otu_map=rules.Filter_otu2seq_map.output
    output: "24.Filter_fasta/pathogens.fasta"
    threads: 1
    log: "logs/Filter_fasta/Filter_fasta.log"
    shell:
        """
        filter_fasta.py \
            -f {input.seqs} \
            -o {output} \
            -m  {input.otu_map} \
            > {log} 2>&1
        """
