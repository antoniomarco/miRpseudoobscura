# R

# libraries
library(DESeq2)

load("results/dpse_DGE.Rdata")

# Differential Gene Expression
#
dps_miR_Lai_myRNA <- dps_miR_pre[,12:14]
cData <-data.frame(sex = c("male","male","female"))
rownames(cData) <-colnames(dps_miR_Lai_myRNA)


# PRECURSORS
ddsMatpreBL <-DESeqDataSetFromMatrix(countData = dps_miR_Lai_myRNA, colData = cData, design = ~ sex)
ddsPreBL <-DESeq(ddsMatpreBL,fitType='local')
resPreBL <-results(ddsPreBL, contrast=c("sex","male","female"))

# save Mohammed results
save(resPreBL, file = "results/dpse_DGE_Lai.Rdata")


# Volcano plot
pdf("plots/dpse_DGE_Lai_volcano.pdf")
hits_over_PreB <-resPreBL[which((resPreBL$padj<0.10) & (resPreBL$log2FoldChange >= log(1.25))),]
hits_under_PreB <-resPreBL[which((resPreBL$padj<0.10) & (resPreBL$log2FoldChange <= log(0.8))),]
plot(resPreBL$log2FoldChange, -log10(resPreBL$padj), xlab = "log2 Fold Change", ylab = "-log10(q)", pch = 16, col = "darkgrey", xlim = c(-max(abs(resPreBL$log2FoldChange),na.rm = TRUE),max(abs(resPreBL$log2FoldChange),na.rm = TRUE)))
abline(h=-log10(0.10), v=c(log(0.8),log(1.25)), lty = 'dashed')
points(hits_over_PreB$log2FoldChange, -log10(hits_over_PreB$padj), pch = 16, col = "darkred")
points(hits_under_PreB$log2FoldChange, -log10(hits_under_PreB$padj), pch = 16, col = "darkblue")
dev.off()

# Smear plot of FDR<0.1
pdf("plots/dpse_DGE_Lai_smear.pdf")
plot(log10(resPreBL$baseMean), resPreBL$log2FoldChange, pch = 19, col = 'white', xlab = "Mean Expression (DESeq2)", ylab = "log2 fold-change", axes = FALSE)
abline(h = 0)
points(log10(resPreBL$baseMean), resPreBL$log2FoldChange, pch = 19, col = 'darkgrey')
values_q_10 <- which((resPreBL$padj<0.10) & (abs(resPreBL$log2FoldChange) >= log(1.25)))
points(log10(resPreBL[values_q_10,]$baseMean), resPreBL[values_q_10,]$log2FoldChange, pch = 19, col = 'black')
axis(1, at=-1:3, tcl= 0.2, labels = c("0.1", "1", "10", "100", "1000"))
axis(2, at = c(-6,-4,-2,0,2,4,6), labels = c(-6,-4,-2,0,2,4,6), las = 2)
dev.off()




pdf("plots/datasets_comparison_l2FC.pdf")
merged_set <- merge(as.data.frame(resPreB), as.data.frame(resPreBL), by = 0)
plot(merged_set$log2FoldChange.x, merged_set$log2FoldChange.y, xlab = "log2 fold-change [GSE179989]", ylab = "log2 fold-change [GSE98013]")
abline(lm(merged_set$log2FoldChange.y ~ merged_set$log2FoldChange.x))
dev.off()

summary(lm(merged_set$log2FoldChange.y ~ merged_set$log2FoldChange.x))
# Multiple R-squared:  0.2933,    Adjusted R-squared:  0.2889 
# F-statistic:  66.4 on 1 and 160 DF,  p-value: 9.962e-14


# Exit
q()


 
