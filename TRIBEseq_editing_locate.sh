#!/bin/sh
#########################################################
#    This BASH script was used to analyse the position
#    of the editing sites in the rice genome using :
#    bedtools(v2.27.1)
#    Please free to use other methods and softwares 
#    as substitution
#########################################################

# ---------Update the following varibles------------
ref_genome="/home/ys/Desktop/japonica_genome/ensemble_plant_genome/ensemble_plant_annotation_file/Oryza_sativa.IRGSP_1.0.dna.toplevel.fa"
ref_anno="/home/ys/Desktop/japonica_genome/ensemble_plant_genome/ensemble_plant_annotation_file/Oryza_sativa.IRGSP_1.0.49.gtf"
ref_repeat_anno="/home/ys/Desktop/japonica_genome/ensemble_plant_genome/ensemble_plant_annotation_file/irgsp1_repeat_unit.gff"

utr5="rice_anno_5utr"
utr3="rice_anno_3utr"
cds="rice_anno_cds"
gene="rice_anno_gene"
exon="rice_anno_exon"
intron="rice_anno_intron"
repeat="rice_anno_repeat"
uncommented="rice_anno_uncommented"

ts=20
# use editing sites under defined coverage threshold. 
# For example, here we use 20 as threshold, editing sites with reads coverage more than 20 were used for next analysis.

# ---------End Update variable------------

########################################################################################################################

# build annotation files for different genomic region

echo "building annotation files for different genomic region"

# build annotation file for different genomic region in 1-based bed format
cat $ref_anno | awk '{FS="\t"}{if($3=="three_prime_utr") print $1"\t"$4+1"\t"$5+1"\t"$3"\t"$7"\t"$9}' > $utr3

cat $ref_anno | awk '{FS="\t"}{if($3=="five_prime_utr") print $1"\t"$4+1"\t"$5+1"\t"$3"\t"$7"\t"$9}' > $utr5

cat $ref_anno | awk '{FS="\t"}{if($3=="CDS") print $1"\t"$4+1"\t"$5+1"\t"$3"\t"$7"\t"$9}' > $cds

cat $ref_anno | awk '{FS="\t"}{if($3=="gene") print $1"\t"$4+1"\t"$5+1"\t"$3"\t"$7"\t"$9}' > $gene

cat $ref_anno | awk '{FS="\t"}{if($3=="exon") print $1"\t"$4+1"\t"$5+1"\t"$3"\t"$7"\t"$9}' > $exon

bedtools subtract -a $gene -b $exon > $intron




echo "============================================================================="
echo
echo using $ts as filtering threshold
echo
echo "============================================================================="

########################################################################################################################

mkdir hices_genomic_pos
prefix=`echo *common*dp$ts`
tissue=${prefix%%_*}

# filter the editing sites by choosing A-to-G mutations with annotation on reference(+) strand and T-to-C mutations with annotation on reverse-complement(-) strand

for pos in $(ls rice_anno_*);
	do bedtools intersect -a *common*dp$ts -b $pos -wb | awk '{FS="\t"}{if($4=="A" && $21=="+"||$4=="T" && $21=="-") print $0}'> $tissue"_usual_hices_dp"$ts"_in_"$pos; 
	sed -i '1 i #r1\tr1\tr1\tr1\tr1\tr1\tr1\tr1\tr2\tr2\tr2\tr2\tr2\tr2\tr2\tr2' $tissue"_usual_hices_dp"$ts"_in_"$pos;
	sed -i '2 i #Chromosome\tStart\tEnd\tREF\tALT\tInformation\tCoverage_Depth\tEditing_Efficiency\tChromosome\tStart\tEnd\tREF\tALT\tInformation\tCoverage_Depth\tEditing_Efficiency\tChromosome\tStart\tEnd\tLocation\tStrand\tGenomic_Information' $tissue"_usual_hices_dp"$ts"_in_"$pos;
	
	echo
	echo "============================================================================="
	echo "put "$tissue"_usual_hices_dp"$ts"_in_"$pos "in folder of hices_genomic_pos"
	echo "============================================================================="
	done

# modified the original repeat annotation file to fit the current editing sites coordinates.
sed 's/chr0//' $ref_repeat_anno | sed 's/chr//' | awk '{FS="\t"}{print $1"\t"$4+1"\t"$5+1"\t"$2"\t"$9}'>$repeat

# next check the editing sites that in repeat and no-repeat genomic regions.
bedtools intersect -a *common*dp$ts -b $repeat -wb > $tissue"_usual_hices_dp"$ts"_in_"$repeat
bedtools subtract -a *common*dp$ts -b $tissue"_usual_hices_dp"$ts"_in_"$repeat > $tissue"_usual_hices_dp"$ts"_in_no_repeat"
sed -i '1 i #r1\tr1\tr1\tr1\tr1\tr1\tr1\tr1\tr2\tr2\tr2\tr2\tr2\tr2\tr2\tr2' $tissue"_usual_hices_dp"$ts"_in_"$repeat
sed -i '2 i #Chromosome\tStart\tEnd\tREF\tALT\tInformation\tCoverage_Depth\tEditing_Efficiency\tChromosome\tStart\tEnd\tREF\tALT\tInformation\tCoverage_Depth\tEditing_Efficiency\tChromosome\tStart\tEnd\tLocation\tGenomic_information' $tissue"_usual_hices_dp"$ts"_in_"$repeat

# check the editing sites that in uncommented genomic regions.
bedtools subtract -a *common*dp$ts -b $tissue"_usual_hices_dp"$ts"_in_"$gene > $tissue"_usual_hices_dp"$ts"_in_"$uncommented
sed -i '1 i #r1\tr1\tr1\tr1\tr1\tr1\tr1\tr1\tr2\tr2\tr2\tr2\tr2\tr2\tr2\tr2' $tissue"_usual_hices_dp"$ts"_in_"$uncommented
sed -i '2 i #Chromosome\tStart\tEnd\tREF\tALT\tInformation\tCoverage_Depth\tEditing_Efficiency\tChromosome\tStart\tEnd\tREF\tALT\tInformation\tCoverage_Depth\tEditing_Efficiency' $tissue"_usual_hices_dp"$ts"_in_"$uncommented

# move the annotation files of different genomic region into /hices_genomic_pos/
mv rice_anno_* hices_genomic_pos 
mv $tissue"_usual_hices_dp"$ts"_in_"* hices_genomic_pos
##############################################################################################################################

# next calculate the unique editing sites in each type of genomic region and unique genomic region that marked by editing sites.

cd hices_genomic_pos
for f in $(ls $tissue"_usual_hices_dp"*rice_anno*);
	do cat  $f |grep -v "^#" | awk '{FS="\t"}{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$20}'| sort | uniq > ${f/usual/uniq};
	echo 
	echo "============================================================================="
	echo "put ${f/usual/uniq} in folder of hices_genomic_pos"
	echo "============================================================================="
	done

# summary the number of high-confident editings in utr3, utr5, cds and intron, repesctively
pos=(3utr 5utr cds intron)
for i in ${pos[@]};
	do for f in $(ls $tissue"_usual"*"hices_dp$ts"_in_rice_anno_$i);
		do echo $i >>"data_summary_of_${f%_$i}_genomic"; cat $f | grep -v "^#" | sort | wc -l >> "data_summary_of_${f%_$i}_genomic"
		done
	done

echo "============================================================================="
echo 
echo "Summary the number of high-confident editings in utr3, utr5, cds and intron into: "
echo "data_summary_of_${f%_$i}_genomic"
echo 
echo "============================================================================="
echo "done!"

#################################################################################################
cd ..
mv TRIBEseq_repeat_count.sh hices_genomic_pos
mv TRIBEseq_genomic_position_in_repeat.sh hices_genomic_pos

echo "============================================================================="
echo 
echo -e "Now, run: \ncd ./hices_genomic_pos && bash TRIBEseq_repeat_count.sh "$tissue"_usual_hices_dp"$ts"_in_"$repeat
echo "to count the different types of repeats marked by high-confident editings"
echo 
echo "============================================================================="

echo all done!
