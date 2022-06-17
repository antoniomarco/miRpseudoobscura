# libraries
library(DESeq2)

# Load D.pse paired DGE data
load(file = "results/dpse_DGE.Rdata")
# Load D.pse unpaired DGE data
load(file = "results/dpse_DGE_Lai.Rdata")
# Load D.mel DGE data
load(file = "results/dmel_DGE.Rdata")


# Filter low count reads
# at least one sample per condition with reads
dsp_dme <- dsp_dme[!(rownames(dsp_dme) %in% names(which(rowSums(cbind(rowSums(dsp_dme[,6:7]==0)>0, rowSums(dsp_dme[,8:9]==0)>0))==2))),]
# baseMean of 1 or greater
dsp_dme <- dsp_dme[which(rownames(dsp_dme) %in% rownames(resPreB[resPreB$baseMean>=1,])),]

# Comparing Dmel and Dpse fold changes
pdf("plots/dps_dme_orthologs.pdf")
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
# Residual standard error: 1.403 on 55 degrees of freedom
# Multiple R-squared:  0.7245,    Adjusted R-squared:  0.7195 
# F-statistic: 144.6 on 1 and 55 DF,  p-value: < 2.2e-16
summary(In_neoX_LM)
# Multiple R-squared:  0.1903,    Adjusted R-squared:  0.1427 
# F-statistic: 3.996 on 1 and 17 DF,  p-value: 0.06184



## UNPAIRED SET resPreBL
# reoad fresh dataset
# Load D.pse paired DGE data
load(file = "results/dpse_DGE.Rdata")
# Load D.pse unpaired DGE data
load(file = "results/dpse_DGE_Lai.Rdata")
# Load D.mel DGE data
load(file = "results/dmel_DGE.Rdata")

# FILTER LOW READ COUNTS
# at least one sample per condition with reads (female has only one sample)
dsp_dme <- dsp_dme[!(rownames(dsp_dme) %in% names(which(rowSums(cbind(rowSums(dsp_dme[,12:13]==0)>0, dsp_dme[,14]==0))==2))),]
# baseMean of 1 or greater
dsp_dme <- dsp_dme[rownames(resPreBL[resPreBL$baseMean>=1,]),]


DGE_L_merged <- merge(data.frame(resPreBL), dsp_dme, by = 0)

pdf("plots/dps_dme_orthologs_Lai.pdf")
plot(DGE_L_merged$log2FoldChange,DGE_L_merged$DME_L2FC, col = 'white', axes = FALSE, xlab = "D. pseudoobscura (log2 fold-change)", ylab = "D. melanogaster (log2 fold-change)")
abline(a = 0, b = 1, lty = 'dashed')
dps_in_neoX <- which(DGE_L_merged$Muller %in% c("XR_D"))
dps_not_in_neoX <- which(!(DGE_L_merged$Muller %in% c("XR_D")))
points(DGE_L_merged$log2FoldChange,DGE_L_merged$DME_L2FC, col = 'darkgrey', pch = 19)
points(DGE_L_merged[dps_in_neoX,]$log2FoldChange,DGE_L_merged[dps_in_neoX,]$DME_L2FC, col = 'black', pch = 19)
Not_in_neoX_LM <- lm(DGE_L_merged[dps_not_in_neoX,]$DME_L2FC ~ DGE_L_merged[dps_not_in_neoX,]$log2FoldChange)
In_neoX_LM <- lm(DGE_L_merged[dps_in_neoX,]$DME_L2FC ~ DGE_L_merged[dps_in_neoX,]$log2FoldChange)
abline(Not_in_neoX_LM, col = 'darkgrey')
abline(In_neoX_LM, col = 'black')
axis(1, at=-8:4, tcl= 0.2, labels = -8:4)
axis(2, at = -8:4, labels = -8:4, las = 2)
dev.off()

# Regression tests:
summary(Not_in_neoX_LM)
# Multiple R-squared:  0.334,     Adjusted R-squared:  0.3263 
# F-statistic: 43.13 on 1 and 86 DF,  p-value: 3.727e-09
summary(In_neoX_LM)
# Multiple R-squared:  3.46e-05,  Adjusted R-squared:  -0.04542 
# F-statistic: 0.0007612 on 1 and 22 DF,  p-value: 0.9782





# Supplementary Table Information

# reoad fresh dataset
# Load D.pse paired DGE data
load(file = "results/dpse_DGE.Rdata")
# Load D.pse unpaired DGE data
load(file = "results/dpse_DGE_Lai.Rdata")
# Load D.mel DGE data
load(file = "results/dmel_DGE.Rdata")

# Format
paired.set <- data.frame(resPreB[rownames(dps_miR_pre),c(1,2,6)])
unpaired.set <- data.frame(resPreBL[rownames(dps_miR_pre),c(1,2,6)])
both.sets <- merge(paired.set, unpaired.set, by = 0)
rownames(both.sets) <- both.sets$Row.names
out.table <- merge(dps_miR_pre, both.sets, by = 0)
rownames(out.table) <- out.table$Row.names
out.table <- out.table[,c(1,2,3,4,5,16:18,20:25)]
# Add Dmel folds
resDme <- data.frame(resDme)
out.table.2 <- merge(out.table, resDme, by.x = 'OrthDme', by.y = 0, all.x = TRUE, all.y = FALSE)


# Save
write.table(out.table.2, file = "results/supplementary_table.tab", quote = FALSE, sep = "\t")



q()
