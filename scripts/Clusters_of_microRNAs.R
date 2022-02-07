# libraries
library(DESeq2)

# Load D.pse DEG data
load(file = "results/dpse_DGE.Rdata")
# Load D.mel DEG data
load(file = "results/dmel_DGE.Rdata")


# Generate input for clustered microRNAs
dps_sorted <- dps_miR_pre[order(dps_miR_pre$Chr, dps_miR_pre$Start, dps_miR_pre$End),]
clusters <- rep(0,nrow(dps_sorted))
clust_num <- 1
# Find clusters
for(i in 2:nrow(dps_sorted) ) {
	if(((dps_sorted[i,'Start'] - dps_sorted[i-1,'End']) < 10000) & (dps_sorted[i,'Chr'] == dps_sorted[i-1,'Chr'])) {
		clusters[i-1] <- clust_num
		clusters[i] <- clust_num
	} else {
		clust_num <- clust_num + 1
	}
}
dps_sorted$clusters <- clusters
clusters_id <- unique(clusters[clusters > 0])

dominant_cluster <- vector()
for (j in 1:length(clusters_id)) {
	temp_clust <- resPreB[rownames(dps_sorted[dps_sorted$clusters == clusters_id[j],]),]
	dominant_cluster <- c(dominant_cluster,rownames(temp_clust[temp_clust$baseMean ==  max(temp_clust$baseMean),]))
}
singles <- rownames(dps_sorted[dps_sorted$clusters == 0,])
dps_miR_clus <- dps_miR_pre[c(singles,dominant_cluster),]



# Differential Gene Expression
#
dps_miR_clus_myRNA <- dps_miR_clus[,6:9]
cData <-data.frame(sex = c("female","female","male","male"), batch = c("A","B","A","B"))
rownames(cData) <-colnames(dps_miR_clus_myRNA)

# PRECURSORS
ddsMatpreC <-DESeqDataSetFromMatrix(countData = dps_miR_clus_myRNA, colData = cData, design = ~ batch + sex) # batch control
ddsPreC <-DESeq(ddsMatpreC,fitType='local')
resPreC <-results(ddsPreC, contrast=c("sex","male","female"))

# Reduce number of tests, lower q-values
resPreC[which((resPreC$padj < 0.1) & (abs(resPreC$log2FoldChange)>=-log2(0.8))),]
# 13 SBEmiRs: 8 female, 5 male



# Comparative analysis
#
# Dmel and Dpse ortholog pairs
dsp_dme <- dps_miR_clus[grepl("^dme",dps_miR_clus$OrthDme),]
dsp_dme$DME_L2FC <- resDme[dsp_dme$OrthDme,]$log2FoldChange
dsp_dme$DPS_L2FC <- resPreB[rownames(dsp_dme),]$log2FoldChange
#
# Comparing Dmel and Dpse fold changes
pdf("plots/dps_dme_orthologs_clus.pdf")
plot(dsp_dme$DPS_L2FC,dsp_dme$DME_L2FC, col = 'white', axes = FALSE, xlab = "D. pseudoobscura (log2 fold-change)", ylab = "D. melanogaster (log2 fold-change)")
abline(a = 0, b = 1, lty = 'dashed')
dps_in_neoX <- which(dsp_dme$Muller %in% c("XR_D"))
dps_not_in_neoX <- which(!(dsp_dme$Muller %in% c("XR_D")))
points(dsp_dme$DPS_L2FC,dsp_dme$DME_L2FC, col = 'darkgrey', pch = 19)
points(dsp_dme[dps_in_neoX,]$DPS_L2FC,dsp_dme[dps_in_neoX,]$DME_L2FC, col = 'black', pch = 19)
Not_in_neoX_LM <- lm(dsp_dme[dps_not_in_neoX,]$DME_L2FC ~ dsp_dme[dps_not_in_neoX,]$DPS_L2FC)
In_neoX_LM <- lm(dsp_dme[dps_in_neoX,]$DME_L2FC ~ dsp_dme[dps_in_neoX,]$DPS_L2FC)
abline(Not_in_neoX_LM, col = 'darkgrey')
abline(In_neoX_LM, col = 'black')
axis(1, at=-8:4, tcl= 0.2, labels = -8:4)
axis(2, at = -8:4, labels = -8:4, las = 2)
dev.off()

# Regression tests:
summary(Not_in_neoX_LM)
#  p-value: < 8.442e-11; adjusted R^2 = 0.5872
summary(In_neoX_LM)
# p = p-value: 0.05896; adjusted R^2 = 0.1254


#
# Novel versus conserved D. pseudoobscura
# ONLY TESTIS AND OVARIES
#
# Remove microRNAs not expressed in either testis or ovary
dps_miR_clus_gonads <- dps_miR_clus[(dps_miR_clus$testis_t > 0) | (dps_miR_clus$ovary_t > 0),]







rpm_testis <- dps_miR_clus_gonads$testis_t/(sum(dps_miR_clus_gonads$testis_t)/1e6)
rpm_ovary <- dps_miR_clus_gonads$ovary_t/(sum(dps_miR_clus_gonads$ovary_t)/1e6)
dps_miR_clus_gonads$testis_ovary_ratio <- log2((rpm_testis+1)/(rpm_ovary+1))
# Novel versus conserved D. pseudoobscura
dpse_specific_X <- rownames(dps_miR_clus_gonads[(dps_miR_clus_gonads$DpsSpecific == "yes") & (dps_miR_clus_gonads$Muller %in% c("XL_A","XR_D")),])
dpse_specific_A <- rownames(dps_miR_clus_gonads[(dps_miR_clus_gonads$DpsSpecific == "yes") & !(dps_miR_clus_gonads$Muller %in% c("XL_A","XR_D")),])
dpse_nospecific_X <- rownames(dps_miR_clus_gonads[(dps_miR_clus_gonads$DpsSpecific == "no") & (dps_miR_clus_gonads$Muller %in% c("XL_A","XR_D")),])
dpse_nospecific_A <- rownames(dps_miR_clus_gonads[(dps_miR_clus_gonads$DpsSpecific == "no") & !(dps_miR_clus_gonads$Muller %in% c("XL_A","XR_D")),])
# Boxplot
pdf("plots/novel_versus_conserved_SB_gonads_clus.pdf")
boxplot(dps_miR_clus_gonads[dpse_nospecific_A,]$testis_ovary_ratio, dps_miR_clus_gonads[dpse_nospecific_X,]$testis_ovary_ratio, dps_miR_clus_gonads[dpse_specific_A,]$testis_ovary_ratio, dps_miR_clus_gonads[dpse_specific_X,]$testis_ovary_ratio, axes = FALSE, col = 'white') 
axis(1, at = 1:4, label = c("Conserved (A)", "Conserved (X)", "Novel (A)", "Novel (X)"))
abline(h=0, col = "darkgrey")
axis(2, at = c(-10, -5, 0, 5, 10), label = c(-10, -5, 0, 5, 10), las = 2)
points(c(rep(1,nrow(dps_miR_clus_gonads[dpse_nospecific_A,])), rep(2,nrow(dps_miR_clus_gonads[dpse_nospecific_X,])), rep(3,nrow(dps_miR_clus_gonads[dpse_specific_A,])), rep(4,nrow(dps_miR_clus_gonads[dpse_specific_X,]))), c(dps_miR_clus_gonads[dpse_nospecific_A,]$testis_ovary_ratio, dps_miR_clus_gonads[dpse_nospecific_X,]$testis_ovary_ratio, dps_miR_clus_gonads[dpse_specific_A,]$testis_ovary_ratio, dps_miR_clus_gonads[dpse_specific_X,]$testis_ovary_ratio), pch = 16, col = "grey") 
dev.off()
# Compare novel X versus all other miR
wilcox.test(dps_miR_clus_gonads[c(dpse_specific_A, dpse_nospecific_X, dpse_nospecific_A),]$testis_ovary_ratio, dps_miR_clus_gonads[dpse_specific_X,]$testis_ovary_ratio)
# p-value = 0.1409



q()
