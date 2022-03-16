#!/bin/sh
#########################################################
#    This BASH script was used to identify editing using 
#         Fastp(v0.20.1)
#         Cutadapt(v3.2)
#         Bowtie(v1.3.0)
#         Trimmomatics(v0.36);
#         SAMtools(v1.7-1);
#         Bedtools(v2.27.1)
#         pileup2base_no_strand.pl
#    Please free to use other methods and softwares 
#    for calling SNPs as substitution
#########################################################
# ---------Update the following varibles------------
ref_genome="/home/ys/Desktop/japonica_genome/ensemble_plant_genome/ensemble_plant_annotation_file/Oryza_sativa.IRGSP_1.0.dna.toplevel.fa"
ref_anno="/home/ys/Desktop/japonica_genome/ensemble_plant_genome/ensemble_plant_annotation_file/Oryza_sativa.IRGSP_1.0.49.gtf"
bowtie_indices="/home/ys/Desktop/bowtie_index/bt/riceREF_bowtie"
TRIMMOMATIC_JAR="/home/ys/biosoft/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar"
PICARD_JAR="/home/ys/biosoft/picard/picard_2.23.1/picard.jar"
pileup2base_no_strand="~/biosoft/pileup2base/pileup2base-master/pileup2base_no_strand.pl"

wget https://mirbase.org/ftp/CURRENT/genomes/osa.gff3
sed -i 's/Chr//g' osa.gff3
sed -i 's/MIR/mir/g' osa.gff3
grep "miRNA_primary_transcript" osa.gff3 > osa_primary_MIRNA.gff3
cat osa_primary_MIRNA.gff3 | awk -F '["\t";=]' '{print $1"\t"$4"\t"$5"\t"$7"\t"$14;}' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | grep -v "#" | sort -n > osa_primary_mirna.gff3

primary_mirna_anno="~/Desktop/sRNA_anno/new_anno/osa_primary_mirna.gff3"



wget https://mirbase.org/ftp/CURRENT/mature.fa.gz
gunzip  mature.fa.gz
sed -n '/>osa/{N;p}' mature.fa > Mature_file.txt

rice_mature_mirna="/home/ys/Desktop/sRNA_anno/new_anno/Mature_file.txt"



gunzip hairpin.fa.gz 
sed -i 's/>/\n>/g' hairpin.fa 
sed -i 's/loop/loop\n/g' hairpin.fa 
sed -n '/>osa-MIR/{N;p}' hairpin.fa > rice_hairpin.fa
sed -i 's/ /_/g' rice_hairpin.fa
sed -i 's/U/T/g' rice_hairpin.fa
bowtie-build rice_hairpin.fa /home/ys/Desktop/sRNA_anno/new_anno/bowtie_index_rice_hairpin/rice_hairpin_bowtie

bowtie_primary_smallRNA_indices="/home/ys/Desktop/sRNA_anno/new_anno/bowtie_index_rice_hairpin/rice_hairpin_bowtie"

# ---------End Update variable------------

# trim raw reads 
for i in $(ls *.fq.gz); do  fastp -w 12 –q 25 –u 20 –n 5 -A -L -i ${i} -o ./${i}_fpout.gz -h ${i}_fpout.html -j ${i}_fpout.json ; done
for i in $(ls *_fpout.gz); do cutadapt -a AGATCGGAAGA -O 5 -e 0 ${i} -o ./${i}_cut3;done
for i in $(ls *cut3); do java -jar /home/ys/biosoft/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 12 -phred33 $i $i"_trimout.fq.gz" LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:15; done
for i in $(ls *trimout.fq.gz);do cutadapt -a "A{20}" -m 18 -M 34 $i -o ${i%fq.gz}_18-34.fq.gz;done

# map trimmed reads to rice genome or rice primary miRNA
for i in $(ls *trimout._18-34.fq.gz);do bowtie -q -p 12 -n 1 -e 50 -a -m 1 --best --strata -x $bowtie_indices $i -S ./trim_genome_align14/${i%.fq.gz}_genome.sam;cd trim_genome_align14; samtools view -@ 12 -bSh -q 20 ${i%.fq.gz}_genome.sam > ./${i%.fq.gz}_genome.bam; rm ${i%.fq.gz}_genome.sam; samtools sort -@ 10 ${i%.fq.gz}_genome.bam -o ${i%.fq.gz}_genome.sorted.bam; samtools index ${i%.fq.gz}_genome.sorted.bam;cd ..;done
for i in $(ls *trimout._18-34.fq.gz);do bowtie -q -p 12 -n 1 -e 50 -a -m 1 --best --strata -x $bowtie_primary_smallRNA_indices $i ./trim_genome_align14/${i}_hairpin_align.bwt;done

# miRNA read counts
for i in $(ls *.bwt);do cat $i | cut -f3 | sort | uniq -c | sed 's/MIR/mir/g' | awk '{FS=" "}{print $1"\t"$2}'> ${i%%_raw*}_against_hairpin_count.txt;done

# check base coverages, select A-to-G and T-to-C mutations as potential editing sites
for i in $(ls *sorted.bam);do samtools mpileup -Q 0 --skip-indels -f $ref_genome $i > ${i%.bam}_pileup.txt;done
for i in $(ls *r*_pileup.txt);do perl $pileup2base_no_strand $i 0 ./${i%.txt}"_pile2base_no_strand.txt";done
for i in $(ls *no_strand.txt);do cat $i | awk '{FS="\t"}{if ($3=="A" && $7 > 0 || $3=="T" && $6 > 0) print $0}' > ${i%.txt}_editing.txt;done

# filter base with depth and editing efficiencies
depth=5
efficiency=0.02

for i in $(ls *pile2base_no_strand_editing.txt);do cat $i | awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > ${i%.txt}.bed; done
for i in $(ls *.bed);do bedtools intersect -a $i -b $primary_mirna_anno -wb > ${i%.bed}_MIRNA.txt;done
for i in $(ls *_MIRNA.txt);do cat $i | awk '{if ($4=="A" && $12=="+" && $8 > '$depth' && $8/($5+$8)>'$efficiency') print $0"\t"$8/($5+$8)} {if ($4=="T" && $12=="-" && $7 > '$depth' && $7/($6+$7) > '$efficiency') print $0"\t"$7/($6+$7)}' > ${i%.txt}_filterDepthEfficiency.txt;done
for i in $(ls *DepthEfficiency.txt);do bedtools subtract -a $i -b ADAR_r1_*_filterDepthEfficiency.txt ADAR_r2_*_filterDepthEfficiency.txt > ${i%.txt}_against_adar.txt;done
bedtools intersect -a DB1AR_r1*_filterDepthEfficiency_against_adar.txt -b DB1AR_r2*_filterDepthEfficiency_against_adar.txt -wb | awk -F '["\t";=]' '{print $1"\t"$2"\t"$3"\t"$4"\t"$14"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15"\t"$16"\t"$17"\t"$18"\t"$28"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27}' > common_editing_MIRNA.txt
cat common_editing_MIRNA.txt | awk '{print $10}' | sort -n | uniq | while read id; do grep -i "$id" -A1  $rice_mature_mirna >> common_mirna.fa; done


