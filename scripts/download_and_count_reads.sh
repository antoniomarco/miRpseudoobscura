#!/usr/bin/bash

############################
# DROSOPHILA PSEUDOOBSCURA #
############################

# PAPER'S ORIGINAL DATA
# Download expression data from EBI ENA
curl -O "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/073/SRR15115773/SRR15115773.fastq.gz" # male 1
curl -O "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/074/SRR15115774/SRR15115774.fastq.gz" # male 2
curl -O "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/075/SRR15115775/SRR15115775.fastq.gz" # female 1
curl -O "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/076/SRR15115776/SRR15115776.fastq.gz" # female 2

#
# Remove adapters
cutadapt -a TGGAATTCTCGGGT -j 4 -m 18 -M 26 SRR15115773.fastq.gz | gzip > dpse_male_1_t.fastq.gz
cutadapt -a TGGAATTCTCGGGT -j 4 -m 18 -M 26 SRR15115774.fastq.gz | gzip > dpse_male_2_t.fastq.gz
cutadapt -a TGGAATTCTCGGGT -j 4 -m 18 -M 26 SRR15115775.fastq.gz | gzip > dpse_female_1_t.fastq.gz
cutadapt -a TGGAATTCTCGGGT -j 4 -m 18 -M 26 SRR15115776.fastq.gz | gzip > dpse_female_2_t.fastq.gz
# Remove raw read files
rm SRR151157?.fastq.gz

#
# Download additional ovary and testis expression data
# Source: https://www.ebi.ac.uk/ena/browser/view/PRJNA208680
curl -O "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR902/SRR902010/SRR902010.fastq.gz" # ovaries
curl -O "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR902/SRR902011/SRR902011.fastq.gz" # testes
# Remove adapters
cutadapt -a ATCTCGTATGCCGT -j 4 -m 18 -M 26 SRR902010.fastq.gz | gzip > SRR902010_t.fastq.gz
cutadapt -a ATCTCGTATGCCGT -j 4 -m 18 -M 26 SRR902011.fastq.gz | gzip > SRR902011_t.fastq.gz
# Remove raw read files
rm SRR90201?.fastq.gz
#
# Download independent male and female reads from https://pubmed.ncbi.nlm.nih.gov/29233922/
curl -O "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR546/007/SRR5461097/SRR5461097.fastq.gz"
curl -O "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR546/008/SRR5461098/SRR5461098.fastq.gz"
curl -O "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR546/009/SRR5461099/SRR5461099.fastq.gz"
# Remove adapters
cutadapt -a CTGTAGGCACCA -j 4 -m 18 -M 26 SRR5461097.fastq.gz | gzip > SRR5461097_t.fastq.gz
cutadapt -a TGGAATTCTCGG -j 4 -m 18 -M 26 SRR5461098.fastq.gz | gzip > SRR5461098_t.fastq.gz
cutadapt -a TGGAATTCTCGG -j 4 SRR5461099.fastq.gz | gzip > SRR5461099_temp.fastq.gz # two-step process
cutadapt -u 4 -u -4 -j 4 -m 18 -M 26 SRR5461099_temp.fastq.gz | gzip > SRR5461099_t.fastq.gz
rm SRR5461099_temp.fastq.gz
# Remove raw read files
rm *109?.fastq.gz
#
# Map to reference genome
for file in *_t.fastq.gz
do
	name=${file%".fastq.gz"}
	hisat2 -p 4 -x Dpse_104 -U $file -S ${name}.sam
done
#
# Remove FASTQ files
rm *.fastq.gz
#
#  Read counts
# ALL COMBINED FEATURE featureCounts
featureCounts -t miRNA_primary_transcript -g Name -a datasets/dps_104.gff3 -o datasets/Dpse_miRNA_expression.tab dpse_female_1_t.sam dpse_female_2_t.sam dpse_male_1_t.sam dpse_male_2_t.sam SRR902010_t.sam SRR902011_t.sam SRR5461097_t.sam SRR5461098_t.sam SRR5461099_t.sam
# featureCounts -t miRNA -g Name -a datasets/dps_104.gff3 -o datasets/Dpse_miRNA_expression_arms.tab dpse_female_1_t.sam dpse_female_2_t.sam dpse_male_1_t.sam dpse_male_2_t.sam SRR902010_t.sam SRR902011_t.sam SRR5461097_t.sam SRR5461098_t.sam SRR5461099_t.sam # Not used in current anlaysis, but it can be of interest
#
# Remove temporary map files
rm *.sam


############################
# DROSOPHILA MELANOGASTER  #
############################


#
# Download published Drosophila melanogaster small RNA reads
curl -O "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR016/SRR016854/SRR016854.fastq.gz" # male
curl -O "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR018/SRR018039/SRR018039.fastq.gz" # female
curl -O "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR069/SRR069836/SRR069836.fastq.gz" # male
curl -O "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR069/SRR069837/SRR069837.fastq.gz" # female
#
# trim adaptors
for file in *.fastq.gz
do
	name=${file%".fastq.gz"}
	cutadapt -a CTGTAGGCA -j 4 -m 18 -M 26 $file | gzip > ${name}_t.fastq.gz
done
#
# Mapping to reference genome (dm6)
for file in *_t.fastq.gz;
do
	name=${file%".fastq.gz"}
	hisat2 -p 4 -x dm6 -U $file -S ${name}.sam
done
#
# Read counts
featureCounts -t miRNA_primary_transcript -g Name -a datasets/dme_dm6.gff3 -o datasets/Dmel_miRNA_expression.tab *.sam
#
# clean
rm *.sam *.fastq.gz

# Exit
exit 0
 
