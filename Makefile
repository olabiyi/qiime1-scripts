.PHONY: run summarize

run:
	snakemake -pr --rerun-incomplete --keep-going --cores 10 --use-singularity --use-conda

summarize:
	snakemake -pr --rerun-incomplete --keep-going --cores 10 --use-singularity 19.Summarize_table/otu_table_rare_taxa_filtered.biom.summary
