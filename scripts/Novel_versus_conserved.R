# R

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
boxplot(dpse_cons_chrom_l2FG ~ dps_l2fc_groups, axes = FALSE, ylab = "log2 fold-change", xlab = "", outline = FALSE, col = 'white', ylim = c(min(dpse_cons_chrom_l2FG), max(dpse_cons_chrom_l2FG))) 
axis(1, at = 1:4, label = c("Conserved (A)", "Conserved (X)", "Novel (A)", "Novel (X)"))
abline(h=0, col = "darkgrey")
axis(2, at = c(-5, 0, 5), label = c(-5, 0, 5), las = 2)
xes <- ifelse(dps_l2fc_groups == "Conserved (A)", 1, ifelse(dps_l2fc_groups == "Conserved (X)", 2, ifelse(dps_l2fc_groups == "Novel (A)", 3, ifelse(dps_l2fc_groups == "Novel (X)", 4, 0))))
xes[xes==1] <- jitter(xes[xes==1], factor = 5)
xes[xes==2] <- jitter(xes[xes==2], factor = 4)
xes[xes==3] <- jitter(xes[xes==3])
xes[xes==4] <- jitter(xes[xes==4])
points(xes,dpse_cons_chrom_l2FG, pch = 16, col = "grey")
dev.off()


log2FC <- c(dpse_nospecific_A_l2FC, dpse_nospecific_X_l2FC,dpse_specific_A_l2FC,dpse_specific_X_l2FC)
chromosome <- c(rep("A",length(dpse_nospecific_A_l2FC)), rep("X",length(dpse_nospecific_X_l2FC)), rep("A",length(dpse_specific_A_l2FC)), rep("X",length(dpse_specific_X_l2FC)))
conservation <- c(rep("cons",length(c(dpse_nospecific_A_l2FC,dpse_nospecific_X_l2FC))), rep("novel",length(c(dpse_specific_A_l2FC, dpse_specific_X_l2FC))))
log2FC_lm <- lm(rank(log2FC) ~ 1 + chromosome * conservation)
anova(log2FC_lm)
# Analysis of Variance Table
# 
# Response: rank(log2FC)
#              Df Sum Sq Mean Sq F value    Pr(>F)    
# chromosome               1    174   173.7  0.2267 0.63510  
# conservation             1   4387  4386.9  5.7243 0.01874 *
# chromosome:conservation  1    216   216.2  0.2821 0.59657  
# Residuals               93  71271   766.4             


#
# USING LAI foldChanges
#
load("results/dpse_DGE_Lai.Rdata")

# Load D.pse DEG data, fresh copy again
load(file = "results/dpse_DGE.Rdata")



# FILTER LOW READ COUNTS
# at least one sample per condition with reads (female has only one sample)
dps_miR_pre <- dps_miR_pre[!(rownames(dps_miR_pre) %in% names(which(rowSums(cbind(rowSums(dps_miR_pre[,12:13]==0)>0, dps_miR_pre[,14]==0))==2))),]
# baseMean of 1 or greater
dps_miR_pre <- dps_miR_pre[rownames(resPreBL[resPreBL$baseMean>=1,]),]


# Novel versus conserved D. pseudoobscura
dpse_specific_X <- rownames(dps_miR_pre[(dps_miR_pre$DpsSpecific == "yes") & (dps_miR_pre$Muller %in% c("XL_A","XR_D")),])
dpse_specific_A <- rownames(dps_miR_pre[(dps_miR_pre$DpsSpecific == "yes") & !(dps_miR_pre$Muller %in% c("XL_A","XR_D")),])
dpse_nospecific_X <- rownames(dps_miR_pre[(dps_miR_pre$DpsSpecific == "no") & (dps_miR_pre$Muller %in% c("XL_A","XR_D")),])
dpse_nospecific_A <- rownames(dps_miR_pre[(dps_miR_pre$DpsSpecific == "no") & !(dps_miR_pre$Muller %in% c("XL_A","XR_D")),])
#
dpse_specific_A_l2FC <- resPreBL[dpse_specific_A,]$log2FoldChange
dpse_specific_A_l2FC <- dpse_specific_A_l2FC[!is.na(dpse_specific_A_l2FC)]
dpse_specific_X_l2FC <- resPreBL[dpse_specific_X,]$log2FoldChange
dpse_specific_X_l2FC <- dpse_specific_X_l2FC[!is.na(dpse_specific_X_l2FC)]
dpse_nospecific_A_l2FC <- resPreBL[dpse_nospecific_A,]$log2FoldChange
dpse_nospecific_A_l2FC <- dpse_nospecific_A_l2FC[!is.na(dpse_nospecific_A_l2FC)]
dpse_nospecific_X_l2FC <- resPreBL[dpse_nospecific_X,]$log2FoldChange
dpse_nospecific_X_l2FC <- dpse_nospecific_X_l2FC[!is.na(dpse_nospecific_X_l2FC)]
dpse_cons_chrom_l2FG <- c(dpse_specific_A_l2FC, dpse_specific_X_l2FC, dpse_nospecific_A_l2FC, dpse_nospecific_X_l2FC)
# Groups
dps_l2fc_groups <- c(rep("Novel (A)", length(dpse_specific_A_l2FC)), rep("Novel (X)", length(dpse_specific_X_l2FC)), rep("Conserved (A)", length(dpse_nospecific_A_l2FC)), rep("Conserved (X)", length(dpse_nospecific_X_l2FC)))
pdf("plots/novel_versus_conserved_SB_all_Mohammed.pdf")
boxplot(dpse_cons_chrom_l2FG ~ dps_l2fc_groups, axes = FALSE, ylab = "log2 fold-change", xlab = "", outline = FALSE, col = 'white', ylim = c(min(dpse_cons_chrom_l2FG), max(dpse_cons_chrom_l2FG))) 
axis(1, at = 1:4, label = c("Conserved (A)", "Conserved (X)", "Novel (A)", "Novel (X)"))
abline(h=0, col = "darkgrey")
axis(2, at = c(-5, 0, 5), label = c(-5, 0, 5), las = 2)
xes <- ifelse(dps_l2fc_groups == "Conserved (A)", 1, ifelse(dps_l2fc_groups == "Conserved (X)", 2, ifelse(dps_l2fc_groups == "Novel (A)", 3, ifelse(dps_l2fc_groups == "Novel (X)", 4, 0))))
xes[xes==1] <- jitter(xes[xes==1], factor = 5)
xes[xes==2] <- jitter(xes[xes==2], factor = 4)
xes[xes==3] <- jitter(xes[xes==3])
xes[xes==4] <- jitter(xes[xes==4])
points(xes,dpse_cons_chrom_l2FG, pch = 16, col = "grey")
dev.off()


log2FC <- c(dpse_nospecific_A_l2FC, dpse_nospecific_X_l2FC,dpse_specific_A_l2FC,dpse_specific_X_l2FC)
chromosome <- c(rep("A",length(dpse_nospecific_A_l2FC)), rep("X",length(dpse_nospecific_X_l2FC)), rep("A",length(dpse_specific_A_l2FC)), rep("X",length(dpse_specific_X_l2FC)))
conservation <- c(rep("cons",length(c(dpse_nospecific_A_l2FC,dpse_nospecific_X_l2FC))), rep("novel",length(c(dpse_specific_A_l2FC, dpse_specific_X_l2FC))))
log2FC_lm <- lm(rank(log2FC) ~ 1 + chromosome * conservation)
anova(log2FC_lm)
# Analysis of Variance Table
# 
# Response: rank(log2FC)
#               Df Sum Sq Mean Sq F value  Pr(>F)  
# chromosome                1   8753  8753.1  2.4120 0.12192   
# conservation              1  29901 29901.2  8.2396 0.00452 **
# chromosome:conservation   1   8174  8174.2  2.2525 0.13491   
# Residuals               209 758453  3629.0          




#
# ONLY TESTIS AND OVARIES
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

dps_miR_pre_gonads <- dps_miR_pre

# Novel versus conserved D. pseudoobscura
dpse_specific_X <- rownames(dps_miR_pre_gonads[(dps_miR_pre_gonads$DpsSpecific == "yes") & (dps_miR_pre_gonads$Muller %in% c("XL_A","XR_D")),])
dpse_specific_A <- rownames(dps_miR_pre_gonads[(dps_miR_pre_gonads$DpsSpecific == "yes") & !(dps_miR_pre_gonads$Muller %in% c("XL_A","XR_D")),])
dpse_nospecific_X <- rownames(dps_miR_pre_gonads[(dps_miR_pre_gonads$DpsSpecific == "no") & (dps_miR_pre_gonads$Muller %in% c("XL_A","XR_D")),])
dpse_nospecific_A <- rownames(dps_miR_pre_gonads[(dps_miR_pre_gonads$DpsSpecific == "no") & !(dps_miR_pre_gonads$Muller %in% c("XL_A","XR_D")),])
# Boxplot
pdf("plots/novel_versus_conserved_SB_gonads.pdf")
all_l2FG <- c(dps_miR_pre_gonads[dpse_nospecific_A,]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_nospecific_X,]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_specific_A,]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_specific_X,]$testis_ovary_ratio)
boxplot(dps_miR_pre_gonads[dpse_nospecific_A,]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_nospecific_X,]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_specific_A,]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_specific_X,]$testis_ovary_ratio, axes = FALSE, col = 'white', outline = FALSE, ylim = c(min(all_l2FG), max(all_l2FG))) 
axis(1, at = 1:4, label = c("Conserved (A)", "Conserved (X)", "Novel (A)", "Novel (X)"))
abline(h=0, col = "darkgrey")
axis(2, at = c(-10, -5, 0, 5, 10), label = c(-10, -5, 0, 5, 10), las = 2)
points(c(jitter(rep(1,nrow(dps_miR_pre_gonads[dpse_nospecific_A,])), factor = 5), jitter(rep(2,nrow(dps_miR_pre_gonads[dpse_nospecific_X,])), factor = 4), jitter(rep(3,nrow(dps_miR_pre_gonads[dpse_specific_A,]))), jitter(rep(4,nrow(dps_miR_pre_gonads[dpse_specific_X,])))), c(dps_miR_pre_gonads[dpse_nospecific_A,]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_nospecific_X,]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_specific_A,]$testis_ovary_ratio, dps_miR_pre_gonads[dpse_specific_X,]$testis_ovary_ratio), pch = 16, col = "grey") 
dev.off()


dps_miR_pre_gonads <- dps_miR_pre_gonads[c(dpse_nospecific_A, dpse_nospecific_X,dpse_specific_A,dpse_specific_X),]
dps_miR_pre_gonads$chr <- c(rep("A",length(dpse_nospecific_A)), rep("X",length(dpse_nospecific_X)), rep("A",length(dpse_specific_A)), rep("X",length(dpse_specific_X)))
dps_miR_pre_gonads$cons <- c(rep("cons",length(c(dpse_nospecific_A,dpse_nospecific_X))), rep("novel",length(c(dpse_specific_A, dpse_specific_X))))
clusters_lm <- lm(rank(testis_ovary_ratio) ~ 1 + chr * cons, data = dps_miR_pre_gonads)
anova(clusters_lm)
# Response: rank(testis_ovary_ratio)
#            Df Sum Sq Mean Sq F value    Pr(>F)    
# chr         1  13435   13435  6.5386  0.01146 *  
# cons        1  48123   48123 23.4211 2.97e-06 ***
# chr:cons    1   1614    1614  0.7855  0.37674    
# Residuals 165 339021    2055 



q()
