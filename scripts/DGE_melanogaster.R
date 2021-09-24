#R

# libraries
library(DESeq2)

dme_miR <- read.table("datasets/Dmel_miRNA_expression.tab", header = TRUE, row.names = 1, skip = 1)
dme_miR_exp <- dme_miR[,6:9]
colnames(dme_miR_exp) <- c("male_1", "female_1", "male_2", "female_2")

cData <-data.frame(sex = c("male","female","male","female"), batch = c("A","A","B","B"))
rownames(cData) <-colnames(dme_miR_exp)
# BATCH CONTROL (negligible effect, but include)
ddsMatDme <-DESeqDataSetFromMatrix(countData = dme_miR_exp, colData = cData, design = ~ batch + sex)
ddsDme <-DESeq(ddsMatDme,fitType='local')
resDme <-results(ddsDme, contrast=c("sex","male","female"))

# Load D. pseudoobscura DGE results
load(file = "results/dpse_DGE.Rdata")

# Dmel and Dpse ortholog pairs
dsp_dme <- dps_miR_pre[grepl("^dme",dps_miR_pre$OrthDme),]
dsp_dme$DME_L2FC <- resDme[dsp_dme$OrthDme,]$log2FoldChange
dsp_dme$DPS_L2FC <- resPreB[rownames(dsp_dme),]$log2FoldChange

save(resDme, dsp_dme, file = "results/dmel_DGE.Rdata")

# Exit
q()
