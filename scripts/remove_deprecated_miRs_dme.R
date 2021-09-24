# R

ddd <- read.table("temp.gff3")
miR_demoted <- c("-287", "-288", "-280", "-289")
# remove:
to_be_removed <- vector()
for(i in 1:length(miR_demoted)) {
	to_be_removed <- c(to_be_removed,which(grepl(miR_demoted[i], ddd$V9)))
}
ddd_o <- ddd[-to_be_removed,]
write.table(ddd_o,"dme_dm6_filtered.gff3", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

q()
