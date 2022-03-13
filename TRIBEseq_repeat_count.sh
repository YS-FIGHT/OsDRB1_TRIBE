#!/bin/sh
#########################################################
#    This BASH script was used to calculate the number 
#    of editing sites and the number of unique repeat
#    unit that marked by editing, resulting a summary
#    file: data_summary_of_"input file"
#########################################################

file=$1
repeat1=$file"_repeat_counts"
repeat2=$file"_repeat_counts_2"
repeat_data="data_summary_of_"$file

other=`grep -v "^#\|MuDR\|Mite\|MITE\|LTR\|SINE\|LINE\|En-Spm\|TYPEU\|hAT" $file | wc -l `
echo -e "other_repeat\t"$other >>$repeat1
MT=`grep -E "Mite|MITE" $1|wc -l `
echo -e "mite\t"$MT >>$repeat1
LT=`grep "LTR" $1|wc -l` 
echo -e "ltr\t"$LT >>$repeat1
MD=`grep "MuDR" $1|wc -l `
echo -e "mudr\t"$MD >>$repeat1
LN=`grep "LINE" $1|wc -l`
echo -e "line\t"$LN >>$repeat1
SN=`grep "SINE" $1|wc -l `
echo -e "sine\t"$SN >>$repeat1
ES=`grep "En-Spm" $1|wc -l `
echo -e "En_Spm\t"$ES >>$repeat1
TP=`grep "TYPEU" $1|wc -l`
echo -e "TYPEU\t"$TP >>$repeat1
HT=`grep "hAT" $1|wc -l`
echo -e "hAT\t"$HT >>$repeat1

sed -i '1 i #Types_of_repeat_unit\t#Amount_of_editing_sites_in_each_type_of_repeat_unit' $repeat1

cat $file |grep -v "^#" | awk '{FS="\t"}{print $NF}'| sort  > temp
other=`grep -v "^#\|MuDR\|Mite\|MITE\|LTR\|SINE\|LINE\|En-Spm\|TYPEU\|hAT" temp | sort | uniq | wc -l`
echo -e "unique_other_repeat\t"$other >>$repeat2
MT=`grep -E "Mite|MITE" temp | sort | uniq | wc -l`
echo -e "unique_mite\t"$MT >>$repeat2
LT=`grep "LTR" temp | sort | uniq | wc -l` 
echo -e "unique_ltr\t"$LT >>$repeat2
MD=`grep "MuDR" temp | sort | uniq | wc -l`
echo -e "unique_mudr\t"$MD >>$repeat2
LN=`grep "LINE" temp | sort | uniq | wc -l`
echo -e "unique_lINE\t"$LN >>$repeat2
SN=`grep "SINE" temp | sort | uniq | wc -l`
echo -e "unique_sine\t"$SN >>$repeat2
ES=`grep "En-Spm" temp | sort | uniq | wc -l`
echo -e "unique En-Spm\t"$ES >>$repeat2
TP=`grep "TYPEU" temp | sort | uniq | wc -l`
echo -e "unique_TYPEU\t"$TP >>$repeat2
HT=`grep "hAT" temp | sort | uniq | wc -l`
echo -e "unique_hAT\t"$HT >>$repeat2

sed -i '1 i #Types_of_repeat_unit\t#Amount_of_uniq_repeats' $repeat2

echo ================================================================================================================================
echo "putting the repeat data of " $file " into directory of " $repeat_data
echo ================================================================================================================================

paste $repeat1 $repeat2 > $repeat_data

rm $repeat1 $repeat2 temp
###################################################################################################
echo ================================================================================================================================
echo 
echo -e "Now, run :\nbash TRIBEseq_genomic_position_in_repeat.sh"
echo 
echo ================================================================================================================================

echo done!
#done
