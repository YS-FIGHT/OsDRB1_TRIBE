#!/bin/sh
#########################################################
#    This BASH script was used to calculate the number 
#    of editing sites of different genomic region that 
#    also in repeat unit, resulting a summary file:
#    "file"_repeat
#    data_summary_of_"file"_repeat
#########################################################

repeat="rice_anno_repeat"
pos=(3utr 5utr cds intron gene uncommented)
for i in ${pos[@]};
	do for f in $(ls *usual*hices*$i);
		do bedtools intersect -a $f -b $repeat -wb > $f"_repeat";
			bash TRIBEseq_repeat_count.sh $f"_repeat";
			a=`grep -v "^#" $f | wc -l`
			b=`cat $f"_repeat" | wc -l` 
			echo -e "no_repeat\t"$[$a-$b] >> "data_summary_of_"$f"_repeat"
			echo -e "total_"$i"_repeat\t"$b >> "data_summary_of_"$f"_repeat"
		done
	done	

echo -e "Now, run : \ncd .. && bash TRIBEseq_flank_seq.sh"
