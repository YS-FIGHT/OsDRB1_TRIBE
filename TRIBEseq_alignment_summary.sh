#!/bin/sh
#########################################################
#    This BASH script was used to summary the alignment 
#    results of all replicates AFTER run 
#    "TRIBEseq_align_&_call_snp.sh" 
#    with paired-end fastq files of all replicates
#########################################################
# ---------Update the following varibles------------

ref_genome_fai="/home/ys/Desktop/japonica_genome/ensemble_plant_genome/ensemble_plant_annotation_file/Oryza_sativa.IRGSP_1.0.dna.toplevel.fa.fai"
ref_anno="/home/ys/Desktop/japonica_genome/ensemble_plant_genome/ensemble_plant_annotation_file/Oryza_sativa.IRGSP_1.0.49.gtf"

# ---------End Update variable------------

# build rice karyotype file
cat $ref_genome_fai | head -n 14 | awk '{print $1"\t1\t"$2}'|grep -Ev "Mt|Pt"| sort -n | sed "1 i Chr\tStart\tEnd" > ./karyotype 

cd align_summary
for i in $(ls *_align_summary.txt);do grep "Uniquely mapped reads %" $i | sed "s/^/${i%_summary.txt}\t&/g" | sed 's/  //g' | sed 's/ |//g' | sed 's/%$//g'  >> uniq_map_summary.txt; done
echo done!

cd ../expression_counts
tissue=(leaf root)
for i in ${tissue[@]};
	do featureCounts -T 8 -p -t exon -g gene_id -a $ref_anno -o ./$i"_all_feature.txt" `echo ${i}*_drb1_r1*HQ30.sort.bam ${i}*_drb1_r2*HQ30.sort.bam ${i}*_adar_r1*HQ30.sort.bam ${i}*_adar_r2*HQ30.sort.bam ${i}*_drb1adar_r1*HQ30.sort.bam ${i}*_drb1adar_r2*HQ30.sort.bam`;
	done
rm *.bam