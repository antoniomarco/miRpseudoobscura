#!/usr/bin/bash

# Software requirements (version used in brackets):
# HISAT2 (2.2.1)
#   Index for Drosophila melanogaster dm6 in $HISAT2_INDEXES required
# featureCounts (2.0.2)
# cutadapt (1.18)
#

# R (3.6.0)
# R libraries:
#  DESeq2 (1.24.0)

# # Output folders
# # Commented out by default
# mkdir datasets results plots

# # Download and build D. pseudoobscura genome index
# # Commented out by default
# # Latest genome assembly
# # UCI_Dpse_MV25, NCBI Release 104 (2020)
# curl -O "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/870/125/GCF_009870125.1_UCI_Dpse_MV25/GCF_009870125.1_UCI_Dpse_MV25_genomic.fna.gz"
# # Uncompress genome sequence and build hisat2 index
# gunzip < GCF_009870125.1_UCI_Dpse_MV25_genomic.fna.gz > GCF_009870125.1_UCI_Dpse_MV25_genomic.fna
# hisat2-build GCF_009870125.1_UCI_Dpse_MV25_genomic.fna Dpse_104
# 
# # Generate microRNA annotation gff3
# # Commented out by default
# bash scripts/generate_dpse_miRNA_gff3.sh
# bash scripts/generate_dmel_miRNA_gff3.sh
# 
# # Download expression datasets and generate read counts tables
# # Commented out by default
# #
# bash scripts/download_and_count_reads.sh
# #
# # REMOVE index genome files
# # Commented out by default
# rm Dpse_104.*


# Differential Gene Expression Analysis: Drosophila pseudoobscura
Rscript scripts/DGE_pseudoobscura.R

# Differential Gene Expression Analysis: Drosophila melanogaster
Rscript scripts/DGE_melanogaster.R

# Orthologs Dmel/Dpse l2FC
Rscript scripts/Compare_Dpse_Dmel_orthologs.R

# New versus conserved microRNAs
Rscript scripts/Novel_versus_conserved.R

# Comparison with other datasets
Rscript scripts/Compare_ours_and_Lai_dataset.R



# EXIT
exit 0
