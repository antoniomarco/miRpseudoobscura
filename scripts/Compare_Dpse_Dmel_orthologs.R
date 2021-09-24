# libraries
library(DESeq2)

# Load D.pse DEG data
load(file = "results/dpse_DGE.Rdata")
# Load D.mel DEG data
load(file = "results/dmel_DGE.Rdata")

# Comparing Dmel and Dpse fold changes
pdf("plots/dps_dme_ortholods.pdf")
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
#  p-value: < 2.2e-16; adjusted R^2 = 0.63
summary(In_neoX_LM)
# p = p-value: 0.01442; adjusted R^2 = 0.2175

q()
