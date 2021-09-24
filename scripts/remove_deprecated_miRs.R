# R

ddd <- read.table("temp.gff3")
miR_demoted <- c("-210b", "-2504", "-2512", "-2514-1", "-2514-2", "-2515", "-2516", "-2518-1", "-2518-2", "-2518-3", "-2518-4", "-2526", "-2528-1", "-2528-2", "-2530-1", "-2530-2", "-2530-3", "-2530-4", "-2530-5", "-2531", "-2532", "-2533", "-2534-1", "-2534-2", "-2534-3", "-2534-4", "-2534-5", "-2534-6", "-2534-7", "-2570", "-2573-1", "-2573-2", "-2573-3", "-2573-4", "-2576-1", "-2576-2", "-2576-3", "-2576-4", "-2576-5", "-3983")
# remove:
to_be_removed <- vector()
for(i in 1:length(miR_demoted)) {
	to_be_removed <- c(to_be_removed,which(grepl(miR_demoted[i], ddd$V9)))
}
ddd_o <- ddd[-to_be_removed,]
write.table(ddd_o,"dps_104_filtered.gff3", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

q()

