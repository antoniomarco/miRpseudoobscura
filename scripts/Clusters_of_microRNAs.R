# libraries
library(DESeq2)

# Load D.pse DEG data
load(file = "results/dpse_DGE.Rdata")
# Load D.mel DEG data
load(file = "results/dmel_DGE.Rdata")


# FILTER LOW READ COUNTS
# at least one sample per condition with reads
dps_miR_pre <- dps_miR_pre[!(rownames(dps_miR_pre) %in% names(which(rowSums(cbind(rowSums(dps_miR_pre[,6:7]==0)>0, rowSums(dps_miR_pre[,8:9]==0)>0))==2))),]
# baseMean of 1 or greater
dps_miR_pre <- dps_miR_pre[rownames(resPreB[resPreB$baseMean>=1,]),]
# remove NA, from previous filter analysis
dps_miR_pre <- dps_miR_pre[!grepl("NA",rownames(dps_miR_pre)),]


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
# 12 SBEmiRs: 8 female, 6 male



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
# F-statistic: 84.78 on 1 and 38 DF,  p-value: 3.215e-11
summary(In_neoX_LM)
# F-statistic: 1.278 on 1 and 16 DF,  p-value: 0.275



#
# Novel versus conserved D. pseudoobscura
## ONLY TESTIS AND OVARIES
#
# Reload clean the Dpse data
load(file = "results/dpse_DGE.Rdata")
# Remove microRNAs not expressed in either testis or ovary
dps_miR_pre <- dps_miR_pre[(dps_miR_pre$testis_t > 0) | (dps_miR_pre$ovary_t > 0),]
# baseMean with DESeq2
# Code from https://support.bioconductor.org/p/75244/ 
dps_miR_pre_myRNA <- dps_miR_pre[,10:11]
cData <-data.frame(gonad = c("ovary","testis"))
rownames(cData) <-colnames(dps_miR_pre_myRNA)
ddsMatpreBG <-DESeqDataSetFromMatrix(countData = dps_miR_pre_myRNA, colData = cData, design = ~ gonad) # batch control
ddsMatpreBG <- estimateSizeFactors(ddsMatpreBG)
baseMeans <- rowMeans(counts(ddsMatpreBG, normalized=TRUE))
# Remove baseMean < 1
dps_miR_pre <- dps_miR_pre[names(which(baseMeans>=1)),]

# Log2FOLD CHANGE
dps_miR_pre$rpm_testis <- dps_miR_pre$testis_t/(sum(dps_miR_pre$testis_t)/1e6)
dps_miR_pre$rpm_ovary <- dps_miR_pre$ovary_t/(sum(dps_miR_pre$ovary_t)/1e6)
dps_miR_pre$testis_ovary_ratio <- log2((dps_miR_pre$rpm_testis+1)/(dps_miR_pre$rpm_ovary+1))
dps_miR_pre$testis_ovary_mean <- (dps_miR_pre$rpm_testis + dps_miR_pre$rpm_ovary) / 2



dps_sortedC <- dps_miR_pre[order(dps_miR_pre$Chr, dps_miR_pre$Start, dps_miR_pre$End),]

clustersC <- rep(0,nrow(dps_sortedC))
clust_num <- 1
# Find clusters
for(i in 2:nrow(dps_sortedC) ) {
	if(((dps_sortedC[i,'Start'] - dps_sortedC[i-1,'End']) < 10000) & (dps_sortedC[i,'Chr'] == dps_sortedC[i-1,'Chr'])) {
		clustersC[i-1] <- clust_num
		clustersC[i] <- clust_num
	} else {
		clust_num <- clust_num + 1
	}
}
dps_sortedC$clustersC <- clustersC
clusters_idC <- unique(clustersC[clustersC > 0])

dominant_cluster <- vector()
for (j in 1:length(clusters_idC)) {
	temp_clust <- dps_sortedC[dps_sortedC$clustersC == clusters_idC[j],]
	dominant_cluster <- c(dominant_cluster,rownames(temp_clust[temp_clust$testis_ovary_mean ==  max(temp_clust$testis_ovary_mean),]))
}
singlesC <- rownames(dps_sortedC[dps_sortedC$clustersC == 0,])
dps_miR_clusC <- dps_miR_pre[c(singlesC,dominant_cluster),]


# Novel versus conserved D. pseudoobscura
dpse_specific_X <- rownames(dps_miR_clusC[(dps_miR_clusC$DpsSpecific == "yes") & (dps_miR_clusC$Muller %in% c("XL_A","XR_D")),])
dpse_specific_A <- rownames(dps_miR_clusC[(dps_miR_clusC$DpsSpecific == "yes") & !(dps_miR_clusC$Muller %in% c("XL_A","XR_D")),])
dpse_nospecific_X <- rownames(dps_miR_clusC[(dps_miR_clusC$DpsSpecific == "no") & (dps_miR_clusC$Muller %in% c("XL_A","XR_D")),])
dpse_nospecific_A <- rownames(dps_miR_clusC[(dps_miR_clusC$DpsSpecific == "no") & !(dps_miR_clusC$Muller %in% c("XL_A","XR_D")),])



# Boxplot
pdf("plots/novel_versus_conserved_SB_gonads_clus.pdf")
all_log2FC <- c(dps_miR_clusC[dpse_nospecific_A,]$testis_ovary_ratio, dps_miR_clusC[dpse_nospecific_X,]$testis_ovary_ratio, dps_miR_clusC[dpse_specific_A,]$testis_ovary_ratio, dps_miR_clusC[dpse_specific_X,]$testis_ovary_ratio)
boxplot(dps_miR_clusC[dpse_nospecific_A,]$testis_ovary_ratio, dps_miR_clusC[dpse_nospecific_X,]$testis_ovary_ratio, dps_miR_clusC[dpse_specific_A,]$testis_ovary_ratio, dps_miR_clusC[dpse_specific_X,]$testis_ovary_ratio, axes = FALSE, col = 'white', outline = FALSE, ylim = c(min(all_log2FC), max(all_log2FC))) 
axis(1, at = 1:4, label = c("Conserved (A)", "Conserved (X)", "Novel (A)", "Novel (X)"))
abline(h=0, col = "darkgrey")
axis(2, at = c(-10, -5, 0, 5, 10), label = c(-10, -5, 0, 5, 10), las = 2)
points(c(jitter(rep(1,nrow(dps_miR_clusC[dpse_nospecific_A,]))), jitter(rep(2,nrow(dps_miR_clusC[dpse_nospecific_X,]))), jitter(rep(3,nrow(dps_miR_clusC[dpse_specific_A,]))), jitter(rep(4,nrow(dps_miR_clusC[dpse_specific_X,])))), c(dps_miR_clusC[dpse_nospecific_A,]$testis_ovary_ratio, dps_miR_clusC[dpse_nospecific_X,]$testis_ovary_ratio, dps_miR_clusC[dpse_specific_A,]$testis_ovary_ratio, dps_miR_clusC[dpse_specific_X,]$testis_ovary_ratio), pch = 16, col = "grey") 
dev.off()

dps_miR_clusC <- dps_miR_clusC[c(dpse_nospecific_A, dpse_nospecific_X,dpse_specific_A,dpse_specific_X),]
dps_miR_clusC$chr <- c(rep("A",length(dpse_nospecific_A)), rep("X",length(dpse_nospecific_X)), rep("A",length(dpse_specific_A)), rep("X",length(dpse_specific_X)))
dps_miR_clusC$cons <- c(rep("cons",length(c(dpse_nospecific_A,dpse_nospecific_X))), rep("novel",length(c(dpse_specific_A, dpse_specific_X))))
clusters_lm <- lm(rank(testis_ovary_ratio) ~ 1 + chr + cons + chr * cons, data = dps_miR_clusC)
anova(clusters_lm)
# Response: rank(testis_ovary_ratio)
#            Df Sum Sq Mean Sq F value   Pr(>F)   
# chr         1   1704  1704.2  1.5621 0.213935   
# cons        1   8118  8118.3  7.4413 0.007393 **
# chr:cons    1    342   341.6  0.3131 0.576877   
# Residuals 113 123279  1091.0  



q()
