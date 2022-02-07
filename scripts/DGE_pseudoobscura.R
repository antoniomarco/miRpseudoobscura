# R

# libraries
library(DESeq2)

# INPUT
dps_miR_pre <- read.table("datasets/Dpse_miRNA_expression.tab", header = TRUE, row.names = 1)

# Muller elements
# Information from J.J. Emerson paper (plus personnal communication)
dps_miR_pre$Muller <- ifelse(dps_miR_pre$Chr == "NC_046679.1", "2_E", "0")
dps_miR_pre$Muller <- ifelse(dps_miR_pre$Chr == "NC_046680.1", "3_C", dps_miR_pre$Muller)
dps_miR_pre$Muller <- ifelse(dps_miR_pre$Chr == "NC_046681.1", "4_B", dps_miR_pre$Muller)
dps_miR_pre$Muller <- ifelse(dps_miR_pre$Chr == "NC_046683.1" & dps_miR_pre$Start <= 32623868, "XL_A", dps_miR_pre$Muller)
dps_miR_pre$Muller <- ifelse(dps_miR_pre$Chr == "NC_046683.1" & dps_miR_pre$Start >= 32627869, "XR_D", dps_miR_pre$Muller)

colnames(dps_miR_pre) <- gsub("SRR902011","testis",gsub("SRR902010","ovary",gsub(".sam","",gsub("dpse_","",colnames(dps_miR_pre)))))
colnames(dps_miR_pre) <- gsub("SRR5461097","male_1_Lai",gsub("SRR5461098","male_2_Lai",gsub("SRR5461099","female_Lai",colnames(dps_miR_pre))))
dps_miR_pre_sorted <- dps_miR_pre[order(dps_miR_pre$Chr, dps_miR_pre$Strand,dps_miR_pre$Start),]
dps_miR_pre_sorted$cluster <- rownames(dps_miR_pre_sorted)

# Conserved and unique microRNAs
#
# Dpse specific microRNAs (Suppl Table 4, Mohammed... Lai. Genome Res)
dps_specific <- c("dps-mir-2502", "dps-mir-2503", "dps-mir-2505", "dps-mir-2506", "dps-mir-2507a", "dps-mir-2507b", "dps-mir-2508", "dps-mir-2510", "dps-mir-2511-1", "dps-mir-2511-2", "dps-mir-2513a", "dps-mir-2513b", "dps-mir-2517a-1", "dps-mir-2517a-2", "dps-mir-2517a-3", "dps-mir-2519", "dps-mir-2520", "dps-mir-2521", "dps-mir-2522a", "dps-mir-2522b", "dps-mir-2523", "dps-mir-2524-1", "dps-mir-2524-2", "dps-mir-2525", "dps-mir-2527", "dps-mir-2529", "dps-mir-2536", "dps-mir-2537", "dps-mir-2538", "dps-mir-2539", "dps-mir-2540", "dps-mir-2541", "dps-mir-2542-1", "dps-mir-2542-2", "dps-mir-2542-3", "dps-mir-2543a-1", "dps-mir-2543a-2", "dps-mir-2543b-1", "dps-mir-2543b-2", "dps-mir-2544", "dps-mir-2545a-1", "dps-mir-2545a-2", "dps-mir-2545a-3", "dps-mir-2545a-4", "dps-mir-2545a-5", "dps-mir-2545a-6", "dps-mir-2545b", "dps-mir-2546", "dps-mir-2547", "dps-mir-2548", "dps-mir-2549", "dps-mir-2550", "dps-mir-2551", "dps-mir-2552", "dps-mir-2554", "dps-mir-2555", "dps-mir-2556", "dps-mir-2557", "dps-mir-2558-1", "dps-mir-2558-2", "dps-mir-2560", "dps-mir-2561", "dps-mir-2562", "dps-mir-2564", "dps-mir-2565", "dps-mir-2566a-1", "dps-mir-2566a-2", "dps-mir-2566b", "dps-mir-2567a", "dps-mir-2567b", "dps-mir-2567c", "dps-mir-2568a", "dps-mir-2568b", "dps-mir-2569", "dps-mir-2571", "dps-mir-2572", "dps-mir-2574a-1", "dps-mir-2574a-2", "dps-mir-2574b", "dps-mir-2575", "dps-mir-276c", "dps_1056", "dps_109", "dps_1105", "dps_113", "dps_1171", "dps_123", "dps_131", "dps_14", "dps_1409", "dps_141", "dps_1425", "dps_15", "dps_1563", "dps_1595", "dps_1668", "dps_1710", "dps_1723", "dps_18", "dps_187", "dps_22", "dps_222", "dps_23", "dps_2346", "dps_239", "dps_245", "dps_2466", "dps_2525", "dps_2543", "dps_2568", "dps_2613", "dps_27", "dps_2737", "dps_2791", "dps_297", "dps_30", "dps_3339", "dps_3340", "dps_3362", "dps_3414", "dps_3415", "dps_3416", "dps_3417", "dps_3418", "dps_3529", "dps_3602", "dps_3611", "dps_3613", "dps_37", "dps_3712", "dps_3750", "dps_3821", "dps_3830", "dps_3831", "dps_3832", "dps_3833", "dps_3855", "dps_41", "dps_42", "dps_455", "dps_463", "dps_53", "dps_55", "dps_61", "dps_62", "dps_7", "dps_706", "dps_72", "dps_73", "dps_76", "dps_774", "dps_79", "dps_794", "dps_83", "dps_856", "dps_87", "dps_872", "dps_88", "dps_91", "dps_934", "dps_96", "dps_998")
#
# Dpse with Dmel ortholog (Suppl Table 4, Mohammed... Lai. Genome Res)
dps_orth_dme <- data.frame(dps = c("dps-bantam", "dps-let-7", "dps-mir-1", "dps-mir-10", "dps-mir-100", "dps-mir-1000", "dps-mir-1002", "dps-mir-1006", "dps-mir-1010", "dps-mir-11", "dps-mir-12", "dps-mir-124", "dps-mir-125", "dps-mir-133", "dps-mir-137", "dps-mir-13a", "dps-mir-13b-1", "dps-mir-13b-2", "dps-mir-14", "dps-mir-184", "dps-mir-190", "dps-mir-193", "dps-mir-210a", "dps-mir-219", "dps-mir-2509", "dps-mir-252", "dps-mir-2535", "dps-mir-2553", "dps-mir-2559", "dps-mir-2563", "dps-mir-263a", "dps-mir-263b", "dps-mir-274", "dps-mir-275", "dps-mir-276a", "dps-mir-276b", "dps-mir-277", "dps-mir-278", "dps-mir-279", "dps-mir-281-1", "dps-mir-281-2", "dps-mir-282", "dps-mir-283", "dps-mir-284", "dps-mir-285", "dps-mir-286", "dps-mir-2a-1", "dps-mir-2a-2", "dps-mir-2b-1", "dps-mir-2b-2", "dps-mir-2c", "dps-mir-3", "dps-mir-304", "dps-mir-305", "dps-mir-306", "dps-mir-307a", "dps-mir-307b", "dps-mir-308", "dps-mir-309", "dps-mir-311", "dps-mir-314", "dps-mir-315", "dps-mir-316", "dps-mir-317", "dps-mir-318", "dps-mir-31a", "dps-mir-31b", "dps-mir-33", "dps-mir-34", "dps-mir-375", "dps-mir-4", "dps-mir-5", "dps-mir-6-1", "dps-mir-6-2", "dps-mir-6-3", "dps-mir-7", "dps-mir-79", "dps-mir-8", "dps-mir-87", "dps-mir-927", "dps-mir-929", "dps-mir-92a", "dps-mir-92b", "dps-mir-92c", "dps-mir-932", "dps-mir-957", "dps-mir-958", "dps-mir-959", "dps-mir-964", "dps-mir-965", "dps-mir-968", "dps-mir-969", "dps-mir-970", "dps-mir-971", "dps-mir-980", "dps-mir-981", "dps-mir-986", "dps-mir-987", "dps-mir-988", "dps-mir-989", "dps-mir-993", "dps-mir-994", "dps-mir-995", "dps-mir-9a", "dps-mir-9b", "dps-mir-9c", "dps-mir-iab-4", "dps-mir-iab-8", "dps_115", "dps_134", "dps_246", "dps_248", "dps_249", "dps_250", "dps_3", "dps_3420", "dps_4", "dps_98"), dme = c("dme-bantam", "dme-let-7", "dme-mir-1", "dme-mir-10", "dme-mir-100", "dme-mir-1000", "dme-mir-1002", "dme-mir-1006", "dme-mir-1010", "dme-mir-11", "dme-mir-12", "dme-mir-124", "dme-mir-125", "dme-mir-133", "dme-mir-137", "dme-mir-13a", "dme-mir-13b-1", "dme-mir-13b-2", "dme-mir-14", "dme-mir-184", "dme-mir-190", "dme-mir-193", "dme-mir-210", "dme-mir-219", "dme-mir-962", "dme-mir-252", "dme-mir-2535b", "dme-mir-310", "dme-mir-1004", "dme-mir-960", "dme-mir-263a", "dme-mir-263b", "dme-mir-274", "dme-mir-275", "dme-mir-276a", "dme-mir-276b", "dme-mir-277", "dme-mir-278", "dme-mir-279", "dme-mir-281-1", "dme-mir-281-2", "dme-mir-282", "dme-mir-283", "dme-mir-284", "dme-mir-285", "dme-mir-286", "dme-mir-2a-1", "dme-mir-2a-2", "dme-mir-2b-1", "dme-mir-2b-2", "dme-mir-2c", "dme-mir-3", "dme-mir-304", "dme-mir-305", "dme-mir-306", "dme-mir-307a", "dme-mir-307b", "dme-mir-308", "dme-mir-309", "dme-mir-312", "dme-mir-314", "dme-mir-315", "dme-mir-316", "dme-mir-317", "dme-mir-318", "dme-mir-31a", "dme-mir-31b", "dme-mir-33", "dme-mir-34", "dme-mir-375", "dme-mir-4", "dme-mir-5", "dme-mir-6-1", "dme-mir-6-2", "dme-mir-6-3", "dme-mir-7", "dme-mir-79", "dme-mir-8", "dme-mir-87", "dme-mir-927", "dme-mir-929", "dme-mir-92a", "dme-mir-92b", "dme-mir-311", "dme-mir-932", "dme-mir-957", "dme-mir-958", "dme-mir-959", "dme-mir-964", "dme-mir-965", "dme-mir-968", "dme-mir-969", "dme-mir-970", "dme-mir-971", "dme-mir-980", "dme-mir-981", "dme-mir-986", "dme-mir-987", "dme-mir-988", "dme-mir-989", "dme-mir-993", "dme-mir-994", "dme-mir-995", "dme-mir-9a", "dme-mir-9b", "dme-mir-9c", "dme-mir-iab-4", "dme-mir-iab-8", "dme-mir-9381", "dme-mir-999", "dme-mir-1003", "dme-mir-998", "dme-mir-1017", "dme-mir-1014", "dme-mir-996", "dme-mir-961", "dme-mir-10404", "dme-mir-956"))
rownames(dps_orth_dme) <- dps_orth_dme$dps
#
dps_miR_pre$OrthDme <- dps_orth_dme[rownames(dps_miR_pre),]$dme
dps_miR_pre$DpsSpecific <- ifelse(rownames(dps_miR_pre) %in% dps_specific, "yes", "no")



# Differential Gene Expression
#
dps_miR_pre_myRNA <- dps_miR_pre[,6:9]
cData <-data.frame(sex = c("female","female","male","male"), batch = c("A","B","A","B"))
rownames(cData) <-colnames(dps_miR_pre_myRNA)

# PRECURSORS
ddsMatpreB <-DESeqDataSetFromMatrix(countData = dps_miR_pre_myRNA, colData = cData, design = ~ batch + sex) # batch control
ddsPreB <-DESeq(ddsMatpreB,fitType='local')
resPreB <-results(ddsPreB, contrast=c("sex","male","female"))



# Volcano plot
pdf("plots/dpse_DGE_volcano.pdf")
hits_over_PreB <-resPreB[which((resPreB$padj<0.10) & (resPreB$log2FoldChange >= log(1.25))),]
hits_under_PreB <-resPreB[which((resPreB$padj<0.10) & (resPreB$log2FoldChange <= log(0.8))),]
plot(resPreB$log2FoldChange, -log10(resPreB$padj), xlab = "log2 Fold Change", ylab = "-log10(q)", pch = 16, col = "darkgrey", xlim = c(-max(abs(resPreB$log2FoldChange),na.rm = TRUE),max(abs(resPreB$log2FoldChange),na.rm = TRUE)))
abline(h=-log10(0.10), v=c(log(0.8),log(1.25)), lty = 'dashed')
points(hits_over_PreB$log2FoldChange, -log10(hits_over_PreB$padj), pch = 16, col = "darkred")
points(hits_under_PreB$log2FoldChange, -log10(hits_under_PreB$padj), pch = 16, col = "darkblue")
dev.off()

# Smear plot of FDR<0.1
pdf("plots/dpse_DGE_smear.pdf")
plot(log10(resPreB$baseMean), resPreB$log2FoldChange, pch = 19, col = 'white', xlab = "Mean Expression (DESeq2)", ylab = "log2 fold-change", axes = FALSE)
abline(h = 0)
points(log10(resPreB$baseMean), resPreB$log2FoldChange, pch = 19, col = 'darkgrey')
values_q_10 <- which((resPreB$padj<0.10) & (abs(resPreB$log2FoldChange) >= log(1.25)))
points(log10(resPreB[values_q_10,]$baseMean), resPreB[values_q_10,]$log2FoldChange, pch = 19, col = 'black')
axis(1, at=-1:3, tcl= 0.2, labels = c("0.1", "1", "10", "100", "1000"))
axis(2, at = c(-6,-4,-2,0,2,4,6), labels = c(-6,-4,-2,0,2,4,6), las = 2)
dev.off()


# Save results as R objects
save(dps_miR_pre, resPreB, file = "results/dpse_DGE.Rdata")

# END
q()

