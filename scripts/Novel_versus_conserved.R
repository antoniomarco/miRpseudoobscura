# R

library(DESeq2)

# Load D.pse DEG data
load(file = "results/dpse_DGE.Rdata")
# Load D.mel DEG data
load(file = "results/dmel_DGE.Rdata")

#
# ALL microRNAs
#
# Novel versus conserved D. pseudoobscura
dpse_specific_X <- rownames(dps_miR_pre[(dps_miR_pre$DpsSpecific == "yes") & (dps_miR_pre$Muller %in% c("XL_A","XR_D")),])
dpse_specific_A <- rownames(dps_miR_pre[(dps_miR_pre$DpsSpecific == "yes") & !(dps_miR_pre$Muller %in% c("XL_A","XR_D")),])
dpse_nospecific_X <- rownames(dps_miR_pre[(dps_miR_pre$DpsSpecific == "no") & (dps_miR_pre$Muller %in% c("XL_A","XR_D")),])
dpse_nospecific_A <- rownames(dps_miR_pre[(dps_miR_pre$DpsSpecific == "no") & !(dps_miR_pre$Muller %in% c("XL_A","XR_D")),])
#
dpse_specific_A_l2FC <- resPreB[dpse_specific_A,]$log2FoldChange
dpse_specific_A_l2FC <- dpse_specific_A_l2FC[!is.na(dpse_specific_A_l2FC)]
dpse_specific_X_l2FC <- resPreB[dpse_specific_X,]$log2FoldChange
dpse_specific_X_l2FC <- dpse_specific_X_l2FC[!is.na(dpse_specific_X_l2FC)]
dpse_nospecific_A_l2FC <- resPreB[dpse_nospecific_A,]$log2FoldChange
dpse_nospecific_A_l2FC <- dpse_nospecific_A_l2FC[!is.na(dpse_nospecific_A_l2FC)]
dpse_nospecific_X_l2FC <- resPreB[dpse_nospecific_X,]$log2FoldChange
dpse_nospecific_X_l2FC <- dpse_nospecific_X_l2FC[!is.na(dpse_nospecific_X_l2FC)]
dpse_cons_chrom_l2FG <- c(dpse_specific_A_l2FC, dpse_specific_X_l2FC, dpse_nospecific_A_l2FC, dpse_nospecific_X_l2FC)
# Groups
dps_l2fc_groups <- c(rep("Novel (A)", length(dpse_specific_A_l2FC)), rep("Novel (X)", length(dpse_specific_X_l2FC)), rep("Conserved (A)", length(dpse_nospecific_A_l2FC)), rep("Conserved (X)", length(dpse_nospecific_X_l2FC)))
pdf("plots/novel_versus_conserved_SB_all.pdf")
boxplot(dpse_cons_chrom_l2FG ~ dps_l2fc_groups, axes = FALSE, ylab = "log2 fold-change", xlab = "") 
axis(1, at = 1:4, label = c("Conserved (A)", "Conserved (X)", "Novel (A)", "Novel (X)"))
abline(h=0, col = "darkgrey")
axis(2, at = c(-5, 0, 5), label = c(-5, 0, 5), las = 2)
points(ifelse(dps_l2fc_groups == "Conserved (A)", 1, ifelse(dps_l2fc_groups == "Conserved (X)", 2, ifelse(dps_l2fc_groups == "Novel (A)", 3, ifelse(dps_l2fc_groups == "Novel (X)", 4, 0)))),dpse_cons_chrom_l2FG, pch = 16, col = "grey")
dev.off()
#
# ONLY TESTIS AND OVARIES
#
# Remove microRNAs not expressed in either testis or ovary
dps_miR_pre_gonads <- dps_miR_pre[(dps_miR_pre$testis_t > 0) | (dps_miR_pre$ovary_t > 0),]
rpm_testis <- dps_miR_pre_gonads$testis_t/(sum(dps_miR_pre_gonads$testis_t)/1e6)
rpm_ovary <- dps_miR_pre_gonads$ovary_t/(sum(dps_miR_pre_gonads$ovary_t)/1e6)
dps_miR_pre_gonads$testis_ovary_ratio <- log2((rpm_testis+1)/(rpm_ovary+1))
# Novel versus conserved D. pseudoobscura
dpse_specific_X <- rownames(dps_miR_pre_gonads[(dps_miR_pre_gonads$DpsSpecific == "yes") & (dps_miR_pre_gonads$Muller %in% c("XL_A","XR_D")),])
dpse_specific_A <- rownames(dps_miR_pre_gonads[(dps_miR_pre_gonads$DpsSpecific == "yes") & !(dps_miR_pre_gonads$Muller %in% c("XL_A","XR_D")),])
dpse_nospecific_X <- rownames(dps_miR_pre_gonads[(dps_miR_pre_gonads$DpsSpecific == "no") & (dps_miR_pre_gonads$Muller %in% c("XL_A","XR_D")),])
dpse_nospecific_A <- rownames(dps_miR_pre_gonads[(dps_miR_pre_gonads$DpsSpecific == "no") & !(dps_miR_pre_gonads$Muller %in% c("XL_A","XR_D")),])
# Boxplot
pdf("plots/novel_versus_conserved_SB_gonads.pdf")
boxplot(dps_miR_pre_gonads[dpse_nospecific_A,]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_nospecific_X,]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_specific_A,]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_specific_X,]$testis_ovary_ratio, axes = FALSE, col = 'white') 
axis(1, at = 1:4, label = c("Conserved (A)", "Conserved (X)", "Novel (A)", "Novel (X)"))
abline(h=0, col = "darkgrey")
axis(2, at = c(-10, -5, 0, 5, 10), label = c(-10, -5, 0, 5, 10), las = 2)
points(c(rep(1,nrow(dps_miR_pre_gonads[dpse_nospecific_A,])), rep(2,nrow(dps_miR_pre_gonads[dpse_nospecific_X,])), rep(3,nrow(dps_miR_pre_gonads[dpse_specific_A,])), rep(4,nrow(dps_miR_pre_gonads[dpse_specific_X,]))), c(dps_miR_pre_gonads[dpse_nospecific_A,]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_nospecific_X,]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_specific_A,]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_specific_X,]$testis_ovary_ratio), pch = 16, col = "grey") 
dev.off()
# Compare novel X versus all other miR
wilcox.test(dps_miR_pre_gonads[c(dpse_specific_A, dpse_nospecific_X, dpse_nospecific_A),]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_specific_X,]$testis_ovary_ratio)
# p-value = 3.781e-05



q()
