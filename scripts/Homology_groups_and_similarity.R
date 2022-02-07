## Load libraries
# install.packages("devtools")
# library(devtools)
# devtools::install_github("mhahsler/rBLAST")
library(rBLAST) # local BLAST wrapper
library(igraph)
library(ips) # Alignment and phylogenetic software wrapper
# devtools::install_github("PhanstielLab/bedtoolsr")
library(bedtoolsr) # Bedtools wrapper

## BLAST database Dmel miRNA precursors
makeblastdb("datasets/dme_dm6.fas", dbtype = "nucl", args = "-out dmel_dm6")
dmel.dm6.db <- blast(db="dmel_dm6", type = "blastn")
# Input FASTA sequences
dps.104.fa <- readDNAStringSet("datasets/dps_104.fas")
dme.dm6.fa <- readDNAStringSet("datasets/dme_dm6.fas")
dme.dps.fa <- c(dps.104.fa,dme.dm6.fa)
# BLAST
output.blast <- predict(dmel.dm6.db, dps.104.fa, BLAST_args = "-evalue 0.1 -word_size 10")

# GRAPH, extract graphs from forrest
blast.graph <- graph_from_data_frame(output.blast[,c(1,2)], directed = FALSE)
blast.components <- components(blast.graph)
more.than.2.homologs <- as.vector(which(table(as.vector(blast.components$membership))>2))

# NOTE RAW to avoid multiple correction errros, but may neet to revise this
for ( i in 1:length(more.than.2.homologs) ) {
	output.alignment <- mafft(as.DNAbin(dme.dps.fa[names(which(blast.components$membership == more.than.2.homologs[i]))]))
	out.image = paste("plots/homology_cluster_", i, ".pdf", sep = "")
	pdf(out.image)
	par(mfrow = c(2,1))
	image.DNAbin(output.alignment, show.bases = FALSE, cex.lab=0.7)
	plot(nj(dist.dna(output.alignment, model = "RAW", pairwise.deletion = FALSE)), "u", , cex = 0.6)
	dev.off()
}
# Only two clusters contained Dpse microrNAs unmatched to putative Dmel orthologs: dps-mir-92c and dps-mir-287. But Dmel has a orthoglogous mir-287 not picked up after filtering
# Dmel mir-92a/b in 3R (E)
# Dpse mir-92a/b in 2 (E) mir-92c in 3 (C). See Maria Ninova paper on the origin of the 311 cluster.
# https://rnajournal.cshlp.org/content/20/3/360.long

## BLAST database Dpse genome
makeblastdb("GCF_009870125.1_UCI_Dpse_MV25_genomic.fna", dbtype = "nucl", args = "-out dps_104")
dps_104.db <- blast(db="dps_104", type = "blastn")
# Input FASTA sequences
# BLAST
output.genome.blast <- predict(dps_104.db, dme.dm6.fa, BLAST_args = "-evalue 0.01 -word_size 10")
# At least 50 nucleotides in the alignment
in.bed <- output.genome.blast[output.genome.blast$Alignment.Length>=50,c(2,9,10,1)]
in.bed$score <- rep(".", nrow(in.bed))
in.bed$strand <- ifelse(in.bed$S.start<in.bed$S.end, "+", "-")
in.bed[,2:3] <- cbind(as.vector(apply(in.bed[,2:3],1,min)),as.vector(apply(in.bed[,2:3],1,max)))
in.bed <- in.bed[order(in.bed$SubjectID, in.bed$S.start, in.bed$S.end),]
number.of.hits <- data.frame(table(in.bed$QueryID))
#
#
# One HIT, DMEL, not annotated in DPSE
#
dmel.onehit.bed <- in.bed[in.bed$QueryID %in% as.vector(number.of.hits[number.of.hits$Freq == 1, 1]),]
novel.dpse <- bedtoolsr::bt.intersect(dmel.onehit.bed,dps.104.precursors.gff3, v = TRUE)
# Potentially newly found orthologs (1 to 1). Out of the scope to evaluate one by one, but we can see the synthesny:
dme.dm6.gff3 <- read.table("datasets/dme_dm6.gff3")
rownames(dme.dm6.gff3) <- gsub("Name=", "", dme.dm6.gff3$V9)
novel.dpse$chr_dme <- dme.dm6.gff3[novel.dpse$V4 ,1]
# Muller elements
novel.dpse$chr_dps <- rep(0, nrow(novel.dpse))
novel.dpse$chr_dps <- ifelse(novel.dpse$V1 == "NC_046679.1", "2_E", "0")
novel.dpse$chr_dps <- ifelse(novel.dpse$V1 == "NC_046680.1", "3_C", novel.dpse$chr_dps)
novel.dpse$chr_dps <- ifelse(novel.dpse$V1 == "NC_046681.1", "4_B", novel.dpse$chr_dps)
novel.dpse$chr_dps <- ifelse(novel.dpse$V1 == "NC_046683.1" & novel.dpse$V2 <= 32623868, "XL_A", novel.dpse$chr_dps)
novel.dpse$chr_dps <- ifelse(novel.dpse$V1 == "NC_046683.1" & novel.dpse$V2 >= 32627869, "XR_D", novel.dpse$chr_dps)
rownames(novel.dpse) <- novel.dpse$V4
novel.dpse <- novel.dpse[,c(1,2,3,6,7,8)]
colnames(novel.dpse) <- c("Scaffold", "Start", "End", "Strand", "Dmel CHR", "Dpse CHR")
write.table(novel.dpse, file = "results/dmel_miR_blast_single_hits_dpse_genome.tab")
# Two incongruent!
# 23 NC_046680.1 13009749 13009804  dme-mir-4945  .  +   chr3R     3_C # In any case, this is between Autosomes. No demasculinization.
# 43 NC_046683.1 37674259 37674347       dme_424  .  +    chrX    XR_D # Both intronic to gene NFAT. Near centromere, either annotation issue or centromeric translocation in Dpse. In any case, X to X.
#
#
# TWO HITS, potential for copies (retroposed or not)
#
two.hits <- number.of.hits[number.of.hits$Freq == 2, ]
two.hits.out <- in.bed[in.bed$QueryID %in% two.hits$Var1,]
# Not already described in Dpse
dps.104.gff3 <- read.table("datasets/dps_104.gff3")
dps.104.precursors.gff3 <- dps.104.gff3[dps.104.gff3$V3 == 'miRNA_primary_transcript', c(1,4,5,9)]
dps.104.precursors.gff3$V9 <- gsub("Name=", "", dps.104.precursors.gff3$V9)
dps.104.precursors.gff3$score <- rep(".",nrow(dps.104.precursors.gff3))
dps.104.precursors.gff3$strand <- ifelse(dps.104.precursors.gff3$V4<dps.104.precursors.gff3$V5, "+", "-")
dps.104.precursors.gff3[,2:3] <- cbind(as.vector(apply(dps.104.precursors.gff3[,2:3],1,min)),as.vector(apply(dps.104.precursors.gff3[,2:3],1,max)))
dps.104.precursors.gff3 <- dps.104.precursors.gff3[order(dps.104.precursors.gff3$V1, dps.104.precursors.gff3$V4, dps.104.precursors.gff3$V5),]
two.hits.out <- bedtoolsr::bt.intersect(two.hits.out,dps.104.precursors.gff3, v = TRUE)
two.hits.out <- two.hits.out[order(two.hits.out$V4),]
# Muller elements
two.hits.out$chr_dme <- dme.dm6.gff3[two.hits.out$V4 ,1]
two.hits.out$chr_dps <- rep(0, nrow(two.hits.out))
two.hits.out$chr_dps <- ifelse(two.hits.out$V1 == "NC_046679.1", "2_E", "0")
two.hits.out$chr_dps <- ifelse(two.hits.out$V1 == "NC_046680.1", "3_C", two.hits.out$chr_dps)
two.hits.out$chr_dps <- ifelse(two.hits.out$V1 == "NC_046681.1", "4_B", two.hits.out$chr_dps)
two.hits.out$chr_dps <- ifelse(two.hits.out$V1 == "NC_046683.1" & two.hits.out$V2 <= 32623868, "XL_A", two.hits.out$chr_dps)
two.hits.out$chr_dps <- ifelse(two.hits.out$V1 == "NC_046683.1" & two.hits.out$V2 >= 32627869, "XR_D", two.hits.out$chr_dps)
two.hits.out <- two.hits.out[,c(4,1,2,3,6,7,8)]
colnames(two.hits.out) <- c("Query", "Scaffold", "Start", "End", "Strand", "Dmel CHR", "Dpse CHR")
write.table(two.hits.out, file = "results/dmel_miR_blast_two_hits_dpse_genome.tab")
# 4 cases. For two the hits are in the same Muller elelemt. For mir-1006 the two hits are in autosomes. For dme_86 the hits are in the Muller element (putative ortholog) and in the XL. If confirmed, that will be in any case a copy form A to X (no demasculinization!)
#
#
# OVER TWO HITS
#
over.two.hits <- number.of.hits[number.of.hits$Freq > 2, ]
colnames(over.two.hits) <- c("Query", "Hits")
write.table(over.two.hits, file = "results/dmel_miR_blast_over_two_hits_dpse_genome.tab")
# rRNA associated microRNAs. Multiple copies in Dmel and Dpse (55). See https://pubmed.ncbi.nlm.nih.gov/25605965/
# dme_208, dme_417, dme_455, dme_5, dme_7, dme-mir-10404
# Probable transposon
# 64     dme-mir-2281   5
# 21 hits in pseudoobscura, but many more in shorter sequences
# Many hits in persimils, in unplaced scaffolds
# The Dmel sequence is also probably NOT a microRNA
# Low complexity
# 121    dme-mir-4973   4
#

# EXIT
q()
