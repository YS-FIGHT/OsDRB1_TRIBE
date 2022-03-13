#!/bin/sh
#########################################################################
# This script is to identify potential and high-confident editing sites using:
#         bcftools(v1.7);
#         bedtools(v2.27.1);
#    Please free to use other methods and softwares 
#    for calling SNPs as substitution.
#########################################################################

#########################################################################

# ---------Update the following varibles------------

input_RBAD1=$1
input_RBAD2=$2
input_RB1=$3
input_RB2=$4
input_AD1=$5
input_AD2=$6
tissue=${1%%_*}

mv $input_RBAD1 $tissue"_drb1adar_against_ref_r1.vcf.gz"
mv $input_RBAD2 $tissue"_drb1adar_against_ref_r2.vcf.gz"
mv $input_RB1 $tissue"_drb1_against_ref_r1.vcf.gz"
mv $input_RB2 $tissue"_drb1_against_ref_r2.vcf.gz"
mv $input_AD1 $tissue"_adar_against_ref_r1.vcf.gz"
mv $input_AD2 $tissue"_adar_against_ref_r2.vcf.gz"

RBAD1=$tissue"_drb1adar_against_ref_r1.vcf.gz"
RBAD2=$tissue"_drb1adar_against_ref_r2.vcf.gz"
RB1=$tissue"_drb1_against_ref_r1.vcf.gz"
RB2=$tissue"_drb1_against_ref_r2.vcf.gz"
AD1=$tissue"_adar_against_ref_r1.vcf.gz"
AD2=$tissue"_adar_against_ref_r2.vcf.gz"

prefix=${RBAD1%_against_ref_r1.vcf.gz}"_against_ref_drb1_adar"
DB1AR_r1=$prefix"_r1"
DB1AR_r2=$prefix"_r2"

# ---------End Update variable------------

# use vcf files generated from previous step
# remove SNPs that exist in both exprimental group(OsDRB1-ADARdd-OE) and control group(ADARdd-OE & OsDRB1-OE), each genotype using 2 replicates
for i in $(ls *.vcf.gz);
	do tabix -p vcf $i;bcftools stats $i > $i.stats
	done
mkdir $DB1AR_r1
mkdir $DB1AR_r2

bcftools isec -p ./$DB1AR_r1/ -C $RBAD1 $RB1 $RB2 $AD1 $AD2 && cd $DB1AR_r1 && mv 0000.vcf ../${DB1AR_r1}".vcf" && cd ../

bcftools isec -p ./$DB1AR_r2/ -C $RBAD2 $RB1 $RB2 $AD1 $AD2 && cd $DB1AR_r2 && mv 0000.vcf ../${DB1AR_r2}".vcf" && cd ../
rm -r $DB1AR_r1 $DB1AR_r2

echo
echo ========================================================================
echo "remove background SNPs that in control lines from DRB1-ADARdd-OE lines"
echo ========================================================================
echo

for i in $(ls *ref_drb1_adar*.vcf);
	do bcftools stats $i > $i.stats
	done
# summarize the nucleotide mutation information
for f in $(ls *stats); 
do grep "^ST" $f | awk '{print $3"\t"$4}'| sed "1 i type_of_mutation\t${f%.vcf.*stats}" > $f"_temp";
done
paste *_temp > $tissue"_replicates_nucleotide_mutation_summary.txt"

# calculating the removed and retained SNPs

paste *drb1adar*ref_drb1_adar_r1*_temp *drb1adar*ref_r1*_temp *drb1adar*ref_drb1_adar_r2*_temp *drb1adar*ref_r2*_temp | awk '{print $1"\t"$2"\t"$4-$2"\t"$6"\t"$8-$6}'| sed '1d' | sed '1 i type_of_mutation\tr1_retained\tr1_removed\tr2_retained\tr2_removed' > $prefix"_mutations_retained_vs_removed.txt"

mkdir nucleotide_mutation_summary && mv *mutation*.txt nucleotide_mutation_summary
rm *temp

# extract RNA nucleotide variation from vcf files
for i in $(ls *.vcf);
	do bcftools view -v snps $i > ${i%.vcf*}"_nucleotide_mutations.vcf";
	done

echo
echo ========================================================================
echo "get nucleotide mutations in DRB1-ADARdd-OE lines"
echo ========================================================================
echo

# extract potential editing sites from RNA nucleotide variation files
for file in $(ls *nucleotide_mutations.vcf);
	do bcftools  filter -O v -o ${file%_nucleotide_mutations*}"_potential_editing_dp0.vcf" -i "( REF = 'A' ) & (ALT='G') | (REF='T' ) & (ALT='C')" $file;
	done

echo
echo ========================================================================
echo "get potential editing sites in DRB1-ADARdd-OE lines"
echo ========================================================================
echo

# filter potential editing sites with different reads coverage thresholds
threshold=(10 20 30)

for file in $(ls *potential_editing_dp0.vcf);
	do for ts in ${threshold[@]};
		do bcftools filter -O v -o ${file%0.vcf*}$ts".vcf" -e "DP < $ts" $file;
		done
	done

echo
echo ==================================================================================
echo "filter potential editing sites in DRB1-ADARdd-OE lines with different thresholds"
echo ==================================================================================
echo

# calculate the editing efficiency of each editing site identified by different reads coverage thresholds

for file in $(ls *dp*0.vcf);
	do bcftools query -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\t%INFO\t%DP\t%DP4{3}\t%DP4{2}\t%DP4{0}\t%DP4{1}\t%DP4{2}\t%DP4{3}\n" $file | awk '{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"($8+$9)/($10+$11+$12+$13)}'| awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}'> $file"_efficiency";
	done
for file in $(ls *efficiency);
	do  sed -i '1 i #Chromosome\tStart\tEnd\tREF\tALT\tInformation\tCoverage_Depth\tEditing_Efficiency' $file;
	done

echo
echo ================================================================================================
echo " calculate editing efficiency of editing sites of DRB1-ADARdd-OE lines in different thresholds"
echo ================================================================================================
echo

# check the common editing sites between replicates under different thresholds

threshold=(0 10 20 30)
for ts in ${threshold[@]};
	do bedtools intersect -a *r1*dp$ts*efficiency -b *r2*dp$ts*efficiency -wb >  $tissue"_common_editing_information_dp"$ts;
		sed -i '1 i #r1\tr1\tr1\tr1\tr1\tr1\tr1\tr1\tr2\tr2\tr2\tr2\tr2\tr2\tr2\tr2' $tissue"_common_editing_information_dp"$ts;
		sed -i '2 i #Chromosome\tStart\tEnd\tREF\tALT\tInformation\tCoverage_Depth\tEditing_Efficiency\tChromosome\tStart\tEnd\tREF\tALT\tInformation\tCoverage_Depth\tEditing_Efficiency' $tissue"_common_editing_information_dp"$ts;
		echo ; 
		echo "=================================================================================================="; 
		echo "Common editing information with threshold of "$ts" was summarized in" $tissue"_common_editing_information_dp"$ts;
		echo "==================================================================================================";
		echo ;
	done

# check the replicate-specific editing sites
threshold=(0 10 20 30)
for ts in ${threshold[@]};
	do bedtools subtract -a *r1*dp$ts*efficiency -b *r2*dp$ts*efficiency > $tissue"_r1_specific_editing_information_dp"$ts;
		sed -i '1 i #Chromosome\tStart\tEnd\tREF\tALT\tInformation\tCoverage_Depth\tEditing_Efficiency' $tissue"_r1_specific_editing_information_dp"$ts;
		echo ;
		echo "===============================================================================================================";
		echo $tissue" r1 specific editing information with threshold of "$ts" was summarized in " $tissue"_r1_specific_editing_information_dp"$ts; 
		echo "===============================================================================================================";
		echo ;
		bedtools subtract -a *r2*dp$ts*efficiency -b *r1*dp$ts*efficiency > $tissue"_r2_specific_editing_information_dp"$ts;
		sed -i '1 i #Chromosome\tStart\tEnd\tREF\tALT\tInformation\tCoverage_Depth\tEditing_Efficiency' $tissue"_r2_specific_editing_information_dp"$ts;
		echo ;
		echo "===============================================================================================================";
		echo $tissue" r2 specific editing information with threshold of "$ts" was summarized in " $tissue"_r2_specific_editing_information_dp"$ts;
		echo "===============================================================================================================";
		echo ;
	done
echo
echo ========================================================
echo "check the overlap editing sites between replicates"
echo "check the replicate-specific editing sites"
echo ========================================================
echo
echo "done!"