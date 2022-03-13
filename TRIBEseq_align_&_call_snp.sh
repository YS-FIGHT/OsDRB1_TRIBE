#!/bin/sh
#########################################################
#    This BASH script was used to call SNPs using 
#         Trimmomatics(v0.36);
#         STAR(v2.5.3a);
#         SAMtools(v1.7-1);
#         featureCounts(v2.0.1);
#         Picard(v2.23.1);
#         bcftools(v1.7);
#    Please free to use other methods and softwares 
#    for calling SNPs as substitution
#########################################################

# ---------Update the following varibles------------
ref_genome="/home/ys/Desktop/japonica_genome/ensemble_plant_genome/ensemble_plant_annotation_file/Oryza_sativa.IRGSP_1.0.dna.toplevel.fa"
ref_anno="/home/ys/Desktop/japonica_genome/ensemble_plant_genome/ensemble_plant_annotation_file/Oryza_sativa.IRGSP_1.0.49.gtf"
star_indices="/home/ys/Desktop/japonica_genome/ensemble_plant_genome/rice_ensemble_plant_STAR_index"
TRIMMOMATIC_JAR="/home/ys/biosoft/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar"
PICARD_JAR="/home/ys/biosoft/picard/picard_2.23.1/picard.jar"
# ---------End Update variable------------

# Use the loop below to run the code for multiple fastq files
#filelist=`ls *.fastq`#
#for file in ${filelist[@]} 
#do

######################################################################
#mkdir TRIBEseq_analysis && cd TRIBEseq_analysis
# all raw fastq files for TRIBEseq should in this directory
######################################################################

# the input fastq file should be pair-end sequenced with .fq* prefix, for example r1.fq or r1.fq.gz
file1=$1
file2=$2
fqfile=${file1%_1.fq*}"_TRIBEanalysis"
mkdir $fqfile && mv $file1 $file2 $fqfile && cd $fqfile

prefix=${file1%_1.fq*} 
trim_outfile=$prefix"_trim_out_"?"P" # P for pair-end sequencing

# Trimmomatics is used to remove low quality bases from either end of the reads. Reads with avg quality score of less than 25 are also removed. Please free to use similar software instead of Trimmomatics

java -jar $TRIMMOMATIC_JAR PE -threads 12 -phred33 $1 $2 -baseout $prefix"_trim_out" LEADING:25 TRAILING:25 AVGQUAL:25 

######################################################################

# Align library with STAR

input=$trim_outfile
echo
echo ==============================
echo "using" $input "for alignment"
echo ==============================
echo

STAR  --runThreadN 18 --outFilterMismatchNoverLmax 0.07 --outFileNamePrefix $prefix"_" --outSAMstrandField intronMotif --outFilterMultimapNmax 1  --genomeDir $star_indices --readFilesIn $input
# --outFilterMismatchNoverLmax 0.07: number of mismatches is <= 7% of mapped read length
# --outFilterMultimapNmax 1: output reads that only map to one loci

output=$prefix".sam"
mv $prefix"_Aligned.out.sam" $output

# following code is tested with samtools 1.7.1, you might have to tweak it a bit bases your installed verison of samtools (these flags can be problematic for older version of samtools: -@, -o)
# remove low quality alignment and convert sam to bam

hqbam=$prefix"_HQ30.bam" 
samtools view -@ 8 -bSh -q 30 $output > $hqbam

# sort the bam file before removing duplicates
sort_out=${hqbam%.bam}".sort.bam"
samtools sort -@ 6 $hqbam -o $sort_out
rm $output $hqbam

cd ..
dir=expression_counts
if [ ! -d $dir ]; then
	mkdir $dir
else
	echo ===================================================================================================
	echo $dir " exists and put the gene expression data of " $prefix " into directory of " $dir
	echo ===================================================================================================
fi

# summary the main aligement information
if [ ! -d align_summary ]; then
	mkdir align_summary
else
	echo "============================"
	echo "align_summary exists"
	echo "============================"
fi

cd $fqfile
grep -E "Number of input reads|Average input read length|Uniquely mapped reads number|Uniquely mapped reads %|Average mapped length" $prefix"_Log.final.out" > $prefix"_align_summary.txt"
sed -i "1 i $prefix"_ailgn_infomation"\tresults" $prefix"_align_summary.txt"
mv $prefix"_align_summary.txt" ../align_summary

# sorted bam file is used for assessing gene expression through featurecounts
count_out=$prefix"_featurecount.txt"
samtools index $sort_out
featureCounts -T 8 -p -t exon -g gene_id -a $ref_anno -o ../$dir/$count_out $sort_out

echo
echo =========================================================================
echo "Counting results can be found in $count_out in $dir file"
echo =========================================================================
echo 


# run Picard to remove duplicates
input_for_picard=$sort_out
dupremove_bam=${sort_out%.bam}"_noDup.bam"
java -Xmx4g -jar $PICARD_JAR MarkDuplicates INPUT=$input_for_picard OUTPUT=$dupremove_bam METRICS_FILE=$prefix"dup.txt" VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true TMP_DIR=tmp ASSUME_SORTED=true

# sort the output bam file from picard
nodup_sort_out=$prefix"_HQ30_noDup_sort.bam"
samtools sort -@ 6 $dupremove_bam -o $nodup_sort_out
rm $dupremove_bam

echo
echo ==============================================================
echo "Done with STAR mapping and PCR duplicate removal with PICARD"
echo "created bam file: $sort_out"
echo "created nodup bam file: $nodup_sort_out"
echo "calling SNPs"
echo ==============================================================
echo
######################################################################

mkdir $prefix"_vcf" 
snpfile=$prefix"_against_REF.vcf.gz"

# call SNPs
samtools mpileup -uvf $ref_genome $nodup_sort_out | bcftools call -vm -Oz > ./$prefix"_vcf"/$snpfile

# check the SNPs content
if [ ! -d ../TRIBEseq_vcffile ];then
	mkdir ../TRIBEseq_vcffile
else
	echo "============================"
	echo "TRIBEseq_vcffile exists"
	echo "============================"
fi

cd $prefix"_vcf"
tabix -p vcf $snpfile
ref_stats=${snpfile%.gz}".stats"
bcftools stats $snpfile > $ref_stats

cp $snpfile ../../TRIBEseq_vcffile

# visualization of SNPs
plot-vcfstats $ref_stats -p ./$ref_stats"_plot"
echo
echo ==============================================================
echo "the SNPs in "${prefix}" are stored in "$prefix"_vcf"
echo "done!"
echo ==============================================================
echo

######################################################################
## end of do for fastq file
## done!
