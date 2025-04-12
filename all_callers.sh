#########################################################################################################################
#	This shell script takes input as folder name of Tumor.								#
#	e.g BI_hg38_nonUMI_v1.sh Tumor 									#
#											#
#       PANEL-Clinical pipeline for Somatic, Germline, CNV & MSI								#
#	-Bioinformatics: Date of Implementation - 31 July 2023.							# 
#	Tools Used in this pipeline											#
#	1.  Fastqc													#
#	2.  Trimmomatic													#
#	3.  BWA	mem													#
#	4.  MarkDuplicates												#
#   6.  Samtools													#
#	7.  BQSR													#
#	8.  Mutect2													#
#	9.  Varscan													#
#	10. VarDict													#
#	11. Lofreq													#
#       12. HaplotypeCaller												#
#	13. Freebayes													#
#       14. Platypus													#
#	15. CNVkit													#
#	16. msisensor pro												#
#########################################################################################################################


## Taking user inputs and setting resource paths
input1="$1"     # Taking input of sample folder name as argument
workdir=$(pwd)  # storing path of Working directory
raw_tumor="$workdir/$input1"   # Assigning path of sample
out="$workdir/Output_hg38_nonUMI_${input1}"       # Creating and Assigning path of output folder
db="/mnt/database"              # Assigning database folder path
hg38="$db/GRCH38.P14/hg38.fa"   # Assigning path of human genome reference file

./notify_start.sh "$input1 Run Submitted" "PANEL panel job Sumitted for $input1" 

# Merging fastq files and Creating list file
cd $input1

fn=$(ls | head -1| cut -f 1-3 -d"_")

echo $fn

zcat *R1* > ${fn}_All_R1_001.fastq
zcat *R2* > ${fn}_All_R2_001.fastq

ls *All_R1* > ../list_${input1}

cd $workdir

input2="list_${input1}"

echo "Test list filename" $input2
#Dynamic allocation of cpus

t="$(nproc --all)"      # Fetching total number of cpus
tnp=$((`expr $t - 30`)) # Maximum number of cpus will be utilized
val=$((`expr $tnp / 2`)) # Arithmatic operation on Max cpu count
np=$((`printf "%.*f\n" 0 $val`)) # Number of shared cpu

echo "test" $t $tnp $val $np
# Removing pre-exiting files if any

rm $out/analysis.log
rm $out/parallel_variantCaller

echo "First $input1 Second $input2" # Print inputs as a test

# Creating folders for the intermediate results
mkdir -p $out/$input1/FastQC_output
mkdir -p $out/$input1/trimmomatic_output
mkdir -p $out/$input1/FastQC_output_after-trimmomatic
mkdir -p $out/$input1/alignments_stats
mkdir -p $out/$input1/bqsr
mkdir -p $out/$input1/markdup
mkdir -p $out/variantCaller/Varscan
mkdir -p $out/variantCaller/Vardict
mkdir -p $out/variantCaller/Mutect2
mkdir -p $out/variantCaller/lofreq
mkdir -p $out/variantCaller/HaplotypeCaller
mkdir -p $out/variantCaller/Platypus
mkdir -p $out/variantCaller/Freebayes
mkdir -p $out/$input1/CNV_calling
mkdir -p $out/$input1/MSI_analysis


eval "$(conda shell.bash hook)"         # Setting bash for conda environment
conda activate wgs_gatk4                # Activating conda environment

# Fastqc and Trimmomatics run for Normal
while IFS= read -r line			# While loop to add scripts for each lane in "parallel_fastqc_Normal" 
do					# "parallel_trimmomatic_Normal" and "parallel_after_trimmomatic_Normal"
	echo $line
	header=$(cat $raw_tumor/$line | head -n 1);  # Extracting header from fastq file. Will not work if fastq file
                                                        # is not gunzip file. Change command to cat
	
	id=$(echo $header | cut -f 3-4 -d ":" | sed 's/@//');  # Extracting run ID from fastq file
	
	echo "header content" $header # printing header
	sm=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-3 -d"_");  # Extracting file name from file list
	echo $id
	echo $sm
done < $input2
# QC checking

echo "Fastqc for $input1 Started" >> $out/analysis.log   # Writing in log file
date >> $out/analysis.log                                # Writing in log file

fastqc $raw_tumor/${sm}_All_R1*  -t $tnp -o $out/$input1/FastQC_output/

fastqc $raw_tumor/${sm}_All_R2*  -t $tnp -o $out/$input1/FastQC_output/

echo "Fastqc for $input1 Completed" >> $out/analysis.log         # Writing in log file
date >> $out/analysis.log                                        # Writing in log file
echo "##############################" >> $out/analysis.log       # Writing in log file

# Adaptor Trimming and cleaning

echo "Trimmomatic for $input1 Started" >> $out/analysis.log      # Writing in log file
date >> $out/analysis.log                                        # Writing in log file

trimmomatic PE -threads $tnp -phred33 $raw_tumor/${sm}_All_R1* $raw_tumor/${sm}_All_R2* $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R1-trimmed_UP.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_UP.fastq ILLUMINACLIP:$db/Trimmomatic_adaptors/adaptors:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50 -trimlog $out/$input1/trimmomatic_output/${sm}_trimlog.txt

echo "Trimmomatic for $input1 Completed" >> $out/analysis.log      # Writing in log file
date >> $out/analysis.log                                        # Writing in log file
echo "##############################" >> $out/analysis.log       # Writing in log file

# QC-rechecking after trimming
echo "Fastqc after trimmomatic for $input1 Started" >> $out/analysis.log         # Writing in log file
date >> $out/analysis.log                                        # Writing in log file

fastqc -t $tnp $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq -o $out/$input1/FastQC_output_after-trimmomatic/

echo "Fastqc after trimmomatic for $input1 Completed" >> $out/analysis.log # Writing in log file
date >> $out/analysis.log                                        # Writing in log file
echo "##############################" >> $out/analysis.log       # Writing in log file

#Fastqc_summary

unzip $out/$input1/FastQC_output_after-trimmomatic/${sm}_R1-trimmed_P_fastqc.zip -d $out/$input1/FastQC_output_after-trimmomatic/
unzip $out/$input1/FastQC_output_after-trimmomatic/${sm}_R2-trimmed_P_fastqc.zip -d $out/$input1/FastQC_output_after-trimmomatic/

if
        grep 'PASS	Basic Statistics' $out/$input1/FastQC_output_after-trimmomatic/${sm}_R1-trimmed_P_fastqc/summary.txt && grep 'PASS	Basic Statistics' $out/$input1/FastQC_output_after-trimmomatic/${sm}_R2-trimmed_P_fastqc/summary.txt
then
        echo Fastqc is pass   >> $out/FASTQC_validation.log 
else
        echo Fastqc failed    >> $out/FASTQC_validation.log
        exit 1
fi

# BWA and read group attachment
echo "BWA for $input1 Started" >> $out/analysis.log      # Writing in log file
date >> $out/analysis.log                                        # Writing in log file
echo "##############################" >> $out/analysis.log       # Writing in log file

while IFS= read -r line			# While loop to add scripts for each lane in "parallel_BWA" 
do	

	echo $line
	header=$(cat $raw_tumor/$line | head -n 1);	# Extracting header from fastq file. Will not work if fastq file
                                                        # is not gunzip file. Change command to cat

	id=$(echo $header | cut -f 3-4 -d ":" | sed 's/@//');	# Extracting run ID from fastq file
	echo "header content" $header	# printing header
	sm=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-3 -d"_");	# Extracting file name from file 
	nam=$(echo $raw_tumor/$line | xargs -n 1 basename | cut -f 1-3 -d"_");	# list
	echo $id
	echo $sm
	echo $nam

echo "BWA for $input1 Started" >> $out/analysis.log      # Writing in log file
date >> $out/analysis.log                                # Writing in log file

	if [ -z "$(ls -A  $out/$input1/trimmomatic_output/)" ]		# If statement
	then

	bwa mem -t $tnp -M -R "@RG\tID:$id\tPL:ILLUMINA\tLB:PANEL\tSM:$sm\tPI:200" $hg38 $raw_tumor/${sm}_R1_001.fastq $raw_tumor/${sm}_R2_001.fastq 2> $out/$input1/${sm}.align.stderr | samtools sort -@ $tnp -o $out/$input1/${sm}.sorted.bam					# Will execute when 'trimmomatic_output' folder is empty

	else

	bwa mem -t $tnp -M -R "@RG\tID:$id\tPL:ILLUMINA\tLB:PANEL\tSM:$sm\tPI:200" $hg38 $out/$input1/trimmomatic_output/${sm}_R1-trimmed_P.fastq $out/$input1/trimmomatic_output/${sm}_R2-trimmed_P.fastq 2> $out/$input1/${sm}.align.stderr | samtools sort -@ $tnp -o $out/$input1/${sm}.sorted.bam	# Will execute when condition is false

	fi

done < "$input2"

echo "BWA for $input1 Completed" >> $out/analysis.log      # Writing in log file
date >> $out/analysis.log                                        # Writing in log file
echo "##############################" >> $out/analysis.log       # Writing in log file

##location of sorted bam 
echo "$out/$input1/${sm}.sorted.bam" >> $out/File_location.log

# Markduplicate
echo "Markduplicate for $input1 Started" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log

java -jar /home/genomics/software/picard-tools-1.119/MarkDuplicates.jar I=$out/$input1/${sm}.sorted.bam O=$out/$input1/markdup/${sm}_dedup.sorted.bam M=$out/$input1/markdup/${sm}_metrics.txt 2> $out/$input1/markdup/${sm}.stderr

echo "Markduplicate for $input1 Completed" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log                                        # Writing in log file
echo "##############################" >> $out/analysis.log       # Writing in log file

echo "Samtool Alignment Statistics for $input1 Started" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log                                        # Writing in log file

sambamba flagstat -t $tnp $out/$input1/${sm}.sorted.bam  > $out/$input1/alignments_stats/${sm}_sorted.txt   

sambamba flagstat -t $tnp  $out/$input1/markdup/${sm}_dedup.sorted.bam > $out/$input1/alignments_stats/${sm}_dedup.sorted.txt	# Getting alignment statistics

echo "Samtool Alignment Statistics for $input1 Completed" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log                                        # Writing in log file
echo "##############################" >> $out/analysis.log       # Writing in log file

# Generating ontarget.bam and indexing it
samtools view -b -h -L /home/genomics/bedfile/PANEL/_complete_bed/hg38.bed $out/$input1/markdup/${sm}_dedup.sorted.bam > $out/$input1/markdup/${sm}_ontarget.bam

samtools index -@ $tnp $out/$input1/markdup/${sm}_ontarget.bam

sammamba flagstat $out/$input1/markdup/${sm}_ontarget.bam > $out/$input1/alignments_stats/${sm}_ontarget.txt

samtools view -q 1 -F 3840 -L "/home/genomics/bedfile/PANEL/_complete_bed/hg38.bed" -c $out/$input1/markdup/${sm}_ontarget.bam 
> $out/$input1/alignments_stats/${sm}_ontarget_withbed.txt

#Baserecaliberation
echo "Baserecalibrator for $input1 Started" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log                                        # Writing in log file

java -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp -Xms30g -Xmx30g -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar BaseRecalibrator --input $out/$input1/markdup/${sm}_ontarget.bam  --reference $hg38 --known-sites $db/dbSNP_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz --known-sites $db/dbSNP_hg38/Homo_sapiens_assembly38.known_indels.vcf --output $out/$input1/bqsr/${sm}_recal_before_data.table -L /home/genomics/bedfile/PANEL/_complete_bed/hg38.bed -CPB 50000

java -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp -Xms30g -Xmx30g -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar ApplyBQSR -R $hg38 -I $out/$input1/markdup/${sm}_ontarget.bam -bqsr $out/$input1/bqsr/${sm}_recal_before_data.table -O $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam -L /home/genomics/bedfile/PANEL/_complete_bed/hg38.bed -CPB 50000

echo "Baserecalibrator for $input1 Completed" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log                                        # Writing in log file
echo "##############################" >> $out/analysis.log       # Writing in log file

##location of bqsr bam 
echo "$out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam" >> $out/File_location.log

echo "Samtool indexing for $input1 Started" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log                                        # Writing in log file

samtools index -@ $tnp $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam

echo "Samtool indexing for $input1 Completed" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log                                        # Writing in log file
echo "##############################" >> $out/analysis.log       # Writing in log file

echo "Samtools mpileup for $input1 Started" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log                                        # Writing in log file

samtools mpileup -B -f $hg38 $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam > $out/$input1/bqsr/${sm}_mpileup.bqsr.dedup.sorted.bam

echo "Samtools mpileup for $input1 Completed" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log                                        # Writing in log file
echo "##############################" >> $out/analysis.log       # Writing in log file

echo "Varscan, Mutect2, Vardict and lofreq for $input1 Started" >> $out/analysis.log    # Writing in log file 
date >> $out/analysis.log

echo $sm > sample_list.txt

echo "java  -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp -Xms30g -Xmx30g -jar  /home/genomics/anaconda3/envs/wgs_gatk4/share/varscan-2.3.7-4/VarScan.jar mpileup2cns $out/$input1/bqsr/${sm}_mpileup.bqsr.dedup.sorted.bam --vcf-sample-list sample_list.txt --min-coverage 10 --min-reads2 5 --output-vcf 1 --variants --regions-file /home/genomics/bedfile/PANEL/_complete_bed/hg38_varscan.bed > $out/variantCaller/Varscan/${input1}_varscan.vcf" >> $out/parallel_variantCaller

echo "java -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp  -Xms30g -Xmx30g  -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar Mutect2 -R $hg38 -I $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam -L /home/genomics/bedfile/PANEL/_complete_bed/hg38.interval_list --native-pair-hmm-threads $tnp --germline-resource /mnt/database/GRCH38.P14/gnomAD/af-only-gnomad.hg38.vcf.gz -O $out/variantCaller/Mutect2/${input1}_Mutect2.vcf.gz" >> $out/parallel_variantCaller

echo "java  -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp -Xms30g -Xmx30g -jar  /home/genomics/software/VarDictJava/build/libs/VarDict-1.8.3.jar -r 5 -G $hg38 -f 0.01 -N PANEL -b $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam -c 1 -S 2 -E 3 -g 4 /home/genomics/bedfile/PANEL/_complete_bed/hg38_vardict.bed | /home/genomics/software/VarDict/var2vcf_valid.pl -A -N $sm -d 10 > $out/variantCaller/Vardict/${input1}_vardict.vcf" >> $out/parallel_variantCaller

echo "lofreq call-parallel --pp-threads $tnp  --call-indels -l /home/genomics/bedfile/PANEL/_complete_bed/hg38.bed -s -S /mnt/database/dbSNP_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz -f $hg38 -o $out/variantCaller/lofreq/${input1}_lofreq.vcf $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam" >> $out/parallel_variantCaller

parallel -j 1 < $out/parallel_variantCaller

java -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp  -Xms30g -Xmx30g  -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar FilterMutectCalls -R $hg38 -V $out/variantCaller/Mutect2/${input1}_Mutect2.vcf.gz -O $out/variantCaller/Mutect2/${input1}_filtered_Mutect2.vcf.gz

#Reformatting lofreq vcf
Rscript lofreqReformat.R $out/variantCaller/lofreq/${input1}_lofreq.vcf 

echo "Varscan, Mutect2, Vardict and lofreq for $input1 Completed" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log
echo "##############################" >> $out/analysis.log       # Writing in log file

##location of vcf of somatic variant callers
#Mutect2
echo "$out/variantCaller/Mutect2/${input1}_filtered_Mutect2.vcf.gz" >> $out/File_location.log
#Vardict
echo "$out/variantCaller/Vardict/${input1}_vardict.vcf"  >> $out/File_location.log
#Varscan 
echo "$out/variantCaller/Varscan/${input1}_varscan.vcf" >> $out/File_location.log
#lofreq
echo "$out/variantCaller/lofreq/${input1}_lofreq_rf.vcf" >> $out/File_location.log

echo "HaplotypeCaller, Platypus and Freebayes for $input1 Started" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log

java -XX:+UseParallelGC -XX:ParallelGCThreads=$tnp  -Xms30g -Xmx30g  -jar /home/genomics/anaconda3/envs/wgs_gatk4/share/gatk4-4.2.6.1-1/gatk-package-4.2.6.1-local.jar HaplotypeCaller -R $hg38 -I $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam --dbsnp $db/dbSNP_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz -O $out/variantCaller/HaplotypeCaller/${sm}_HaplotypeCaller.vcf  --native-pair-hmm-threads $tnp -L /home/genomics/bedfile/PANEL/_complete_bed/hg38.interval_list

conda deactivate

eval "$(conda shell.bash hook)"
conda activate rna_preprocessing

platypus callVariants --nCPU $tnp --bamFiles=$out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam --refFile=$hg38  --regions=/home/genomics/bedfile/PANEL/_complete_bed/hg38.bed  --output=$out/variantCaller/Platypus/${sm}_platypus.vcf.gz

conda deactivate

zcat $out/variantCaller/Platypus/${sm}_platypus.vcf.gz > $out/variantCaller/Platypus/${sm}_platypus.vcf

sed -i 's/=TC,/=DP,/g' $out/variantCaller/Platypus/${sm}_platypus.vcf
sed -i 's/=TR,/=AD,/g' $out/variantCaller/Platypus/${sm}_platypus.vcf
sed -i 's/##FORMAT=<ID=NR,/##FORMAT=<ID=DP,/g' $out/variantCaller/Platypus/${sm}_platypus.vcf
sed -i 's/##FORMAT=<ID=NV,/##FORMAT=<ID=AD,/g' $out/variantCaller/Platypus/${sm}_platypus.vcf
sed -i 's/;TC=/;DP=/g' $out/variantCaller/Platypus/${sm}_platypus.vcf
sed -i 's/;TR=/;AD=/g' $out/variantCaller/Platypus/${sm}_platypus.vcf
sed -i 's/GT:GL:GOF:GQ:NR:NV/GT:GL:GOF:GQ:DP:AD/g' $out/variantCaller/Platypus/${sm}_platypus.vcf

freebayes -f $hg38 -t /home/genomics/bedfile/PANEL/_complete_bed/hg38.bed -b $out/$input1/bqsr/${sm}_bqsr.dedup.sorted.bam  > $out/variantCaller/Freebayes/${sm}_freebayes.vcf

echo "HaplotypeCaller, Platypus and Freebayes for $input1 Completed" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log

##location of vcf of germline variant callers
#Haplotypecaller
echo "$out/variantCaller/HaplotypeCaller/${sm}_HaplotypeCaller.vcf" >> $out/File_location.log
#Platypus
echo "$out/variantCaller/Platypus/${sm}_platypus.vcf"  >> $out/File_location.log
#Freebayes
echo "$out/variantCaller/Freebayes/${sm}_freebayes.vcf" >> $out/File_location.log


echo "CNV calling for $input1 Started" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log                                        # Writing in log file
echo "##############################" >> $out/analysis.log       # Writing in log file

# Activating conda env for cnvkit
eval "$(conda shell.bash hook)"
conda activate cnvkit

cnvkit.py batch $out/$input1/${sm}.sorted.bam -r /home/genomics/bedfile/CNV/my_reference.cnn -d $out/$input1/CNV_calling/

cnvkit.py segmetrics $out/$input1/CNV_calling/${sm}.sorted.cnr -s $out/$input1/CNV_calling/${sm}.sorted.cns --ci --t-test -o $out/$input1/CNV_calling/${sm}_segment.cns --drop-low-coverage

awk ' NR==1;{ if($11<0.05) {print} } ' $out/$input1/CNV_calling/${sm}_segment.cns > $out/$input1/CNV_calling/${sm}_segment_sig.cns

cnvkit.py call $out/$input1/CNV_calling/${sm}_segment_sig.cns --filter ci -m threshold -o $out/$input1/CNV_calling/${sm}_segment_call.cns

cnvkit.py export vcf $out/$input1/CNV_calling/${sm}_segment_call.cns -i $sm -o $out/$input1/CNV_calling/${sm}_cnv_call_final.vcf

echo "CNV calling for $input1 Completed" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log                                        # Writing in log file
echo "##############################" >> $out/analysis.log     # Writing in log file

# location of vcf of CNV calling
echo "$out/$input1/CNV_calling/${sm}_cnv_call_final.vcf" >> $out/File_location.log

echo "MSI Analysis for $input1 Started" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log                                        # Writing in log file
echo "##############################" >> $out/analysis.log       # Writing in log file

# Activating conda env for MSI
eval "$(conda shell.bash hook)"
conda activate msi

msisensor-pro pro -d /home/genomics/bedfile/MSI/reference_baseline -t $out/$input1/${sm}.sorted.bam  -o $out/$input1/MSI_analysis/${sm} -e /home/genomics/bedfile/PANEL/_MSI_bed/_MSI_only.bed

echo "MSI Analysis for $input1 Completed" >> $out/analysis.log    # Writing in log file
date >> $out/analysis.log                                        # Writing in log file
echo "##############################" >> $out/analysis.log       # Writing in log file

# location of MSI output
echo "$out/$input1/MSI_analysis/${sm}" >> $out/File_location.log

eval "$(conda shell.bash hook)"         # Setting bash for conda environment
conda activate wgs_gatk4                # Activating conda environment

# print the version of the tools in the log file 
fastqc --version >> $out/version.log
echo "##############################" >> $out/version.log
echo "Trimmomatic"   >> $out/version.log
trimmomatic -version >> $out/version.log
echo "##############################" >> $out/version.log
echo "BWA 0.7.17-r1188" >> $out/version.log
echo "##############################" >> $out/version.log
echo "gatk4-4.2.6.1-1" >> $out/version.log
echo "##############################" >> $out/version.log
echo "sambamba 0.8.2" >> $out/version.log
echo "##############################" >> $out/version.log
echo "samtools 1.6" >> $out/version.log
echo "##############################" >> $out/version.log
echo "VarScan v2.3" >> $out/version.log
echo "##############################" >> $out/version.log
echo "VarDict 1.8" >> $out/version.log
echo "##############################" >> $out/version.log
echo  "LoFreq version 2" >> $out/version.log
echo "##############################" >> $out/version.log
echo  "freebayes v1.3.6" >> $out/version.log
echo "##############################" >> $out/version.log
echo  "platypus 0.8.1" >> $out/version.log
echo "##############################" >> $out/version.log
echo  "CNVkit 0.9.10" >> $out/version.log
echo "##############################" >> $out/version.log
echo  "msisensor-pro v1.2.0" >> $out/version.log
echo "##############################" >> $out/version.log

./notify_complete.sh "$input1 Run Completed" "PANEL panel job Completed for $input1" "$out/analysis.log" "$out/File_location.log"

