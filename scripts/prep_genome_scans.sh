#!/bin/bash

# Download papAnu2 coding sequences (CDS) from Ensembl
wget ftp://ftp.ensembl.org/pub/release-86/fasta/papio_anubis/cds/Papio_anubis.PapAnu2.0.cds.all.fa.gz -O data/papAnu2_cds.fa.gz
gunzip data/papAnu2_cds.fa.gz

# Download papAnu2 gene transfer format (GTF) file from Ensembl
wget ftp://ftp.ensembl.org/pub/release-87/gtf/papio_anubis/Papio_anubis.PapAnu2.0.87.chr.gtf.gz -O data/papAnu2.gtf.gz
gunzip data/papAnu2.gtf.gz

# Download papAnu2 top-level genome from Ensembl
wget ftp://ftp.ensembl.org/pub/release-87/fasta/papio_anubis/dna/Papio_anubis.PapAnu2.0.dna.toplevel.fa.gz -O data/papAnu2.fa.gz
gunzip data/papAnu2.fa.gz

# Index genome
samtools faidx data/papAnu2.fa

cut -f 1-2 data/papAnu2.fa.fai | head -n 20 > data/papAnu2_chr_lengths.txt

# Download Mmul8.0.1 proteome from Uniprot
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000006718_9544.fasta.gz -O data/UP000006718_9544.fasta.gz
gunzip data/UP000006718_9544.fasta.gz

# Download PANTHER mappings for macaque genome
wget ftp://ftp.pantherdb.org/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/PTHR11.0_macacque -O data/PTHR11.0_macacque

# Simplify PANTHER mappings to make it more lightweight
cat data/PTHR11.0_macacque | cut -f 1,3,6,7,9,10 | grep '#' | sort -u > data/panther_pathways_short_macaque.txt

# Build CDS info
scripts/prep_cds_info.R data/papAnu2_cds.fa data/papAnu2.gtf

# Fix headers to just the gene ID: sesbio (https://github.com/sestaton/sesbio) must be installed with the gene_annotation folder in the PATH
clean_multifasta.pl -i data/papAnu2_cds.fa -o data/papAnu2_cds_cleanheaders.fa

# Get open reading frames. HMMER2GO (https://github.com/sestaton/HMMER2GO) should be installed and in the PATH
hmmer2go getorf -i data/papAnu2_cds_cleanheaders.fa -o data/genes_orfs.faa

# Map baboon ORFs to the macaque proteome
phmmer -o data/genes_orfs_uniprot.stdout --tblout data/genes_orfs_uniprot.tblout --domtblout data/genes_orfs_uniprot.domtblout --acc data/genes_orfs.faa data/UP000006718_9544.fasta
cat data/genes_orfs_uniprot.tblout | grep '^tr' | tr -s ' ' | cut -d ' ' -f 1,3,5,19 > data/genes_to_uniprot.txt

# Filter results to build PANTHER annotations for baboon genome
scripts/filter_phmmer_panther_results.R

exit
