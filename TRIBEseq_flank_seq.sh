#!/bin/sh
#########################################################
#    This BASH script was used to analyse the flanking
#    sequence of the editing sites in the rice genome using :
#    bedtools(v2.27.1)
#    seqkit(v0.15.0)
#    MEME(v4.11.2)
#    Please free to use other methods and softwares 
#    as substitution
#########################################################

# ---------Update the following varibles------------
ref_genome="/home/ys/Desktop/japonica_genome/ensemble_plant_genome/ensemble_plant_annotation_file/Oryza_sativa.IRGSP_1.0.dna.toplevel.fa"

ts=10	# use 10 as threshold, please edit as needed
# ---------End Update variable------------


for f in $(ls *common_editing*dp$ts);
	do cat $f | awk '{FS="\t"}{if($4=="A") print$0}'> $f"_AG_origin_pos";		
	done
for f in $(ls *common_editing*dp$ts);
	do cat $f | awk '{FS="\t"}{if($4=="T") print$0}'> $f"_TC_origin_pos";
	done
for f in $(ls *common*$ts*origin_pos);
	do awk '{FS="\t"}{print $1"\t"$2-5"\t"$2+6}' $f > ${f%origin_pos}"_5nt_pos";
	done
for f in $(ls *common*$ts*_5nt_pos);
	do bedtools getfasta -fi $ref_genome -bed $f -fo $f"_flank_5nt.fa";
	done 	
echo 
echo "==================================================================================================================="
echo "using "$ts" as reads coverage threshold to analyse the common high-confident editing sites between replicates"
echo "==================================================================================================================="
seqkit seq *$ts*TC*_flank_5nt.fa -r -p > common_editing_information_dp$ts"_TC_pos_5nt_pos_flank_5nt_RC.fa";

tissue=`echo *common_editing*dp$ts`
prefix=${tissue%common_editing_information_dp$ts}
cat *$ts*AG*_flank_5nt.fa | grep -v '>' > $prefix"hices_5nt_flank.fq"
cat *$ts*TC*_flank_5nt_RC.fa | grep -v '>' >> $prefix"hices_5nt_flank.fq"

for f in $(ls *common*$ts*origin_pos);
	do awk '{FS="\t"}{print $1"\t"$2-250"\t"$2+251}' $f > ${f%origin_pos}"_250nt_pos";
	done
for f in $(ls *common*_250nt_pos);
	do bedtools sort -i $f > $f"_sorted";
	done
for f in $(ls *_250nt_pos_sorted);
	do bedtools merge -i $f > $f"_merge";
	done

for p in $(ls *common*$ts*_merge);
	do bedtools getfasta -fi $ref_genome -bed $p -fo $p".fa";
	done

seqkit seq *$ts*TC*_merge.fa -r -p > common_editing_information_dp$ts"_TC_250nt_pos_sorted_merge_RC.fa";

cat *$ts*AG*_250nt_pos_sorted_merge.fa > $prefix"hices_250nt_flank_merge.fa"
cat *$ts*TC*_pos_sorted_merge_RC.fa >> $prefix"hices_250nt_flank_merge.fa"

mkdir hices_flank_seq
mv *AG* hices_flank_seq 
mv *TC* hices_flank_seq 
mv *hices*flank*.* hices_flank_seq
cd *hices_flank_seq && rm *common*nt*

for f in $(ls *hices*nt*.f*);
	do sed -i 's/T/U/g' $f;
	done
cd ..
echo 
echo "===================================================="
echo "motif search with" $prefix"hices_250nt_flank_merge.fa"
echo "===================================================="

meme hices_flank_seq/$prefix"hices_250nt_flank_merge.fa" -oc ./hices_flank_seq/meme_250nt_merge/ -minw 10 -maxw 12 -dna -maxsize 1000000 -nmotifs 12 -mod zoops

echo done!
