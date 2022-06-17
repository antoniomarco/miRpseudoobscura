#!/usr/bin/bash

# miRBase D_mel hairpin sequences
# NOTE The default hairpin file from miRBase truncates the hairpin sequence
curl -O "https://mirbase.org/ftp/CURRENT/genomes/dme.gff3"
cat dme.gff3 | grep miRNA_primary_transcript | awk '{print $9}' | awk -F';' '{print $1}' | sed 's/ID=//g' | while read accnum;
do
 wget -O - https://www.mirbase.org/cgi-bin/get_seq.pl?acc=$accnum | grep -v -e 'pre'
done  | perl -ne '$_ =~ s/U/T/g; print $_;' | awk '{print $1}' > dme_hairpin.fa


# Map hairpins to reference genome
hisat2 dm6 -f dme_hairpin.fa > dme_dm6_map.sam

# Uniquely mapped microRNA precursors:
cat dme_dm6_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '0' {print $3 "\t.\tmiRNA_primary_transcript\t" $4 "\t" $4+length($10)-1 "\t.\t+\t.\tName=" $1}' > dme_dm6_pos.gff3
cat dme_dm6_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '16' {print $3 "\t.\tmiRNA_primary_transcript\t" $4 "\t" $4+length($10)-1 "\t.\t-\t.\tName=" $1}' > dme_dm6_neg.gff3
cat dme_dm6_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '0' {print $3 "\t.\tmiRNA\t" $4 "\t" $4+length($10)/2 "\t.\t+\t.\tName=" $1 "-5p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' > dme_dm6_pos_arm.gff3
cat dme_dm6_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '0' {print $3 "\t.\tmiRNA\t" $4+length($10)/2+1 "\t" $4+length($10)-1 "\t.\t+\t.\tName=" $1 "-3p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' >> dme_dm6_pos_arm.gff3
cat dme_dm6_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '16' {print $3 "\t.\tmiRNA\t" $4+length($10)/2+1 "\t" $4+length($10)-1 "\t.\t-\t.\tName=" $1 "-5p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' > dme_dm6_neg_arm.gff3
cat dme_dm6_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '16' {print $3 "\t.\tmiRNA\t" $4 "\t" $4+length($10)/2 "\t.\t-\t.\tName=" $1 "-3p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' >> dme_dm6_neg_arm.gff3

# Add additional Dpse microRNAs annotated in https://pubmed.ncbi.nlm.nih.gov/29233922/
wget https://genome.cshlp.org/content/suppl/2017/12/12/gr.226068.117.DC1/Supplemental_12flies_website.zip
unzip Supplemental_12flies_website.zip
for file in 12flies/12FliesAlignments/*.html; do cat $file | grep -A 6 -e '/dme' | perl -ne '$_ =~ s/\<\/SPAN\>//g;$_ =~ s/<SPAN class=\"[a-zA-Z]+\"\>//g; print $_;' | grep -e'dme' -e'<td nowrap>' | perl -ne '$_ =~ s/\<.*\"\>//g; $_ =~ s/\<\/.*\>//g; $_ =~ s/\<td nowrap\>\s+//g;print $_;' | awk '{print $1}'; done | perl -ne 'if($_ =~ /^dme/){print ">".$_}else{$_ =~ s/\-//g; print $_;}' | tac | perl -e '$aa=join("",<>); $aa =~ s/\n\>/\t/g; print $aa;' | awk '{print $2 "\t" $1}' | grep -e '^dme_'| sort | uniq | perl -e '%out = ( ); @lines = split("\n",join("\n",<>)); for(@lines){@row = split("\t",$_); $out{$row[0]}=$row[1];}; for(keys %out){print ">".$_."\n".$out{$_}."\n";}' | grep -e'd' -e'T' > dme_newLai_pre.fa
# FASTA file with even NT numbers to avoid errors in rounding numbers in a later step
cat dme_newLai_pre.fa | perl -ne 'if($_ =~ /^>/){print $_;}else{chomp $_; if(length($_) % 2){chop $_; print "$_\n";}else{print "$_\n";}}' > dme_newLai_pre_even.fa
# Map to reference
hisat2 dm6 -f dme_newLai_pre_even.fa > dme_newLai_pre_dm6_map.sam
# Uniquely mapped precursor miRs
cat dme_newLai_pre_dm6_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '0' {print $3 "\t.\tmiRNA_primary_transcript\t" $4 "\t" $4+length($10)-1 "\t.\t+\t.\tName=" $1}' > dme_dm6_L_pos.gff3
cat dme_newLai_pre_dm6_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '16' {print $3 "\t.\tmiRNA_primary_transcript\t" $4 "\t" $4+length($10)-1 "\t.\t-\t.\tName=" $1}' > dme_dm6_L_neg.gff3
cat dme_newLai_pre_dm6_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '0' {print $3 "\t.\tmiRNA\t" $4 "\t" $4+length($10)/2 "\t.\t+\t.\tName=" $1 "-5p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' > dme_dm6_L_pos_arm.gff3
cat dme_newLai_pre_dm6_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '0' {print $3 "\t.\tmiRNA\t" $4+length($10)/2+1 "\t" $4+length($10)-1 "\t.\t+\t.\tName=" $1 "-3p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' >> dme_dm6_L_pos_arm.gff3
cat dme_newLai_pre_dm6_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '16' {print $3 "\t.\tmiRNA\t" $4+length($10)/2+1 "\t" $4+length($10)-1 "\t.\t-\t.\tName=" $1 "-5p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' > dme_dm6_L_neg_arm.gff3
cat dme_newLai_pre_dm6_map.sam | grep -e NH:i:1 | awk -F'\t' '$2 == '16' {print $3 "\t.\tmiRNA\t" $4 "\t" $4+length($10)/2 "\t.\t-\t.\tName=" $1 "-3p"}' | perl -ne '$_ =~ s/-mir-/-miR-/g; print $_;' >> dme_dm6_L_neg_arm.gff3

# Temporal gff annotation file
cat dme_dm6_neg.gff3 dme_dm6_pos.gff3 dme_dm6_pos_arm.gff3 dme_dm6_neg_arm.gff3 dme_dm6_L_neg.gff3 dme_dm6_L_pos.gff3 dme_dm6_L_pos_arm.gff3 dme_dm6_L_neg_arm.gff3 > temp.gff3

# Remove deprecated microRNAs according to https://pubmed.ncbi.nlm.nih.gov/29233922/
Rscript scripts/remove_deprecated_miRs_dme.R

# FINAL GFF3 FILE
echo -e "##gff-version 3\n##date 2021-6-6\n# Chromosomal coordinates of Drosophila melanogaster microRNAs\n# genome-build:  dm6\n#" > dme_dm6_head.gff3
cat dme_dm6_head.gff3 dme_dm6_filtered.gff3 | sort -k 9 > datasets/dme_dm6.gff3

# Fasta file
cat dme_hairpin.fa dme_newLai_pre.fa | grep -v -e 'dme-mir-287' | grep -v -e 'dme-mir-288' | grep -v -e 'dme-mir-280' | grep -v -e 'dme-mir-289' | grep -A 1 -e'>' | grep -v -e '--' > datasets/dme_dm6.fas

# Remove temp files
rm -r dme_hairpin.fa dme_newLai_pre.fa dme_newLai_pre_even.fa *.gff3 *.sam *12flies* hairpin.fa.gz

# EXIT
exit 0
