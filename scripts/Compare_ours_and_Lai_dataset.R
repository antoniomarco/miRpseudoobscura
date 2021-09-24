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

pdf("plots/datasets_comparison_l2FC.pdf")
plot(resPreB$log2FoldChange,resPreBL$log2FoldChange, xlab = "log2 fold-change [GSE179989]", ylab = "log2 fold-change [GSE98013]")
abline(lm(resPreBL$log2FoldChange ~ resPreB$log2FoldChange))
dev.off()

summary(lm(resPreBL$log2FoldChange ~ resPreB$log2FoldChange))



# Exit
q()


 
