#!/usr/bin/bash

# miRBase D_pse hairpin sequences
# NOTE The default hairpin file from miRBase truncates the hairpin sequence
curl -O "https://mirbase.org/ftp/CURRENT/genomes/dps.gff3"
cat dps.gff3 | grep miRNA_primary_transcript | awk '{print $9}' | awk -F';' '{print $1}' | sed 's/ID=//g' | while read accnum;
do
	wget -O - https://www.mirbase.org/cgi-bin/get_seq.pl?acc=$accnum | grep -v -e 'pre'
done  | perl -ne '$_ =~ s/U/T/g; print $_;' | awk '{print $1}' > dps_hairpin.fa


# Map hairpins to reference genome
hisat2 Dpse_104 -f dps_hairpin.fa > dps_104_map.sam

# Uniquely mapped microRNA precursors:
cat dps_104_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '0' {print $3 "\t.\tmiRNA_primary_transcript\t" $4 "\t" $4+length($10)-1 "\t.\t+\t.\tName=" $1}' > dps_104_pos.gff3
cat dps_104_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '16' {print $3 "\t.\tmiRNA_primary_transcript\t" $4 "\t" $4+length($10)-1 "\t.\t-\t.\tName=" $1}' > dps_104_neg.gff3
cat dps_104_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '0' {print $3 "\t.\tmiRNA\t" $4 "\t" $4+length($10)/2 "\t.\t+\t.\tName=" $1 "-5p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' > dps_104_pos_arm.gff3
cat dps_104_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '0' {print $3 "\t.\tmiRNA\t" $4+length($10)/2+1 "\t" $4+length($10)-1 "\t.\t+\t.\tName=" $1 "-3p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' >> dps_104_pos_arm.gff3
cat dps_104_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '16' {print $3 "\t.\tmiRNA\t" $4+length($10)/2+1 "\t" $4+length($10)-1 "\t.\t-\t.\tName=" $1 "-5p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' > dps_104_neg_arm.gff3
cat dps_104_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '16' {print $3 "\t.\tmiRNA\t" $4 "\t" $4+length($10)/2 "\t.\t-\t.\tName=" $1 "-3p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' >> dps_104_neg_arm.gff3

# Add additional Dpse microRNAs annotated in https://pubmed.ncbi.nlm.nih.gov/29233922/
wget https://genome.cshlp.org/content/suppl/2017/12/12/gr.226068.117.DC1/Supplemental_12flies_website.zip
unzip Supplemental_12flies_website.zip
for file in 12flies/12FliesAlignments/*.html; do cat $file | grep -A 6 -e '/dps' | perl -ne '$_ =~ s/\<\/SPAN\>//g;$_ =~ s/<SPAN class=\"[a-zA-Z]+\"\>//g; print $_;' | grep -e'dps' -e'<td nowrap>' | perl -ne '$_ =~ s/\<.*\"\>//g; $_ =~ s/\<\/.*\>//g; $_ =~ s/\<td nowrap\>\s+//g;print $_;' | awk '{print $1}'; done | perl -ne 'if($_ =~ /^dps/){print ">".$_}else{$_ =~ s/\-//g; print $_;}' | tac | perl -e '$aa=join("",<>); $aa =~ s/\n\>/\t/g; print $aa;' | awk '{print $2 "\t" $1}' | grep -e '^dps_'| sort | uniq | perl -e '%out = ( ); @lines = split("\n",join("\n",<>)); for(@lines){@row = split("\t",$_); $out{$row[0]}=$row[1];}; for(keys %out){print ">".$_."\n".$out{$_}."\n";}' | grep -e'd' -e'T' > dps_newLai_pre.fa
# FASTA file with even NT numbers to avoid errors in rounding numbers in a later step
cat dps_newLai_pre.fa | perl -ne 'if($_ =~ /^>/){print $_;}else{chomp $_; if(length($_) % 2){chop $_; print "$_\n";}else{print "$_\n";}}' > dps_newLai_pre_even.fa
# Map to reference
hisat2 Dpse_104 -f dps_newLai_pre_even.fa > dps_newLai_pre_104_map.sam
# Uniquely mapped precursor miRs
cat dps_newLai_pre_104_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '0' {print $3 "\t.\tmiRNA_primary_transcript\t" $4 "\t" $4+length($10)-1 "\t.\t+\t.\tName=" $1}' > dps_104_L_pos.gff3
cat dps_newLai_pre_104_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '16' {print $3 "\t.\tmiRNA_primary_transcript\t" $4 "\t" $4+length($10)-1 "\t.\t-\t.\tName=" $1}' > dps_104_L_neg.gff3
cat dps_newLai_pre_104_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '0' {print $3 "\t.\tmiRNA\t" $4 "\t" $4+length($10)/2 "\t.\t+\t.\tName=" $1 "-5p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' > dps_104_L_pos_arm.gff3
cat dps_newLai_pre_104_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '0' {print $3 "\t.\tmiRNA\t" $4+length($10)/2+1 "\t" $4+length($10)-1 "\t.\t+\t.\tName=" $1 "-3p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' >> dps_104_L_pos_arm.gff3
cat dps_newLai_pre_104_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '16' {print $3 "\t.\tmiRNA\t" $4+length($10)/2+1 "\t" $4+length($10)-1 "\t.\t-\t.\tName=" $1 "-5p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' > dps_104_L_neg_arm.gff3
cat dps_newLai_pre_104_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '16' {print $3 "\t.\tmiRNA\t" $4 "\t" $4+length($10)/2 "\t.\t-\t.\tName=" $1 "-3p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' >> dps_104_L_neg_arm.gff3

# Temporal gff annotation file
cat dps_104_neg.gff3 dps_104_pos.gff3 dps_104_pos_arm.gff3 dps_104_neg_arm.gff3 dps_104_L_neg.gff3 dps_104_L_pos.gff3 dps_104_L_pos_arm.gff3 dps_104_L_neg_arm.gff3 > temp.gff3

# Remove deprecated microRNAs according to https://pubmed.ncbi.nlm.nih.gov/29233922/
Rscript scripts/remove_deprecated_miRs.R

# FINAL GFF3 FILEecho -e "##gff-version 3\n##date 2021-6-6\n# Chromosomal coordinates of Drosophila pseudoobscura microRNAs\n# genome-build-accession:  GCF_009870125.1_UCI_Dpse_MV25_genomic\n#" > dps_104_head.gff3
cat dps_104_head.gff3 dps_104_filtered.gff3 | sort -k 9 > datasets/dps_104.gff3


# Generate fasta file with precursors
cat datasets/dps_104.gff3 | grep primary_transcript | awk '{print $1 "\t" $4 "\t" $5 "\t" $9 "\t.\t" $7}' | sed 's/Name=//g' > dps_104.bed
bedtools getfasta -name -fi GCF_009870125.1_UCI_Dpse_MV25_genomic.fna -bed dps_104.bed > datasets/dps_104.fas

# Remove temp files
rm -r GCF_009870125.1_UCI_Dpse_MV25_genomic.fna dps_hairpin.fa dps_newLai_pre.fa dps_newLai_pre_even.fa *.gff3 *.sam *12flies* hairpin.fa.gz GCF_009870125.1_UCI_Dpse_MV25_genomic.fna.gz dps_104.bed

# EXIT
exit 0
