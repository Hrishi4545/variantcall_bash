#!/bin/bash;

echo "###################################"
echo "####Last modified on 25-05-2022####"
echo "#Last implementation on 10-06-2022#"
echo "##################################"

echo -e "\n\n\n";
echo "#######################################";
echo "#variables and intialization#";
echo "#######################################";

read1="${1}";
read2="${2}";
ref_fa="${3}";
known_variants="${4}";

#Root Directories

pipeline_root="/mnt/e/Genespectrum/sample_data/";
gatk_root="/mnt/e/Genespectrum/tools/gatk-4.2.2.0/";
fastqc_root="/mnt/e/Genespectrum/tools/FastQC/";
snpeff_root="/mnt/e/Genespectrum/tools/snpEff/";
picard_root="/mnt/e/Genespectrum/tools/";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo "#######################################";
echo "#indexing reference genome#";
echo "#######################################";

bwa index -a bwtsw "${3}";

basename "${3}" .fa;

#create sequence dictionary
java -jar "${gatk_root}"gatk.jar CreateSequenceDictionary R="${3}" O="${3}".dict


date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo "#######################################";
echo "#mapping/alignment of reads with reference genome#";
echo "#######################################";

bwa mem -M -t 4 "${3}" "${1}" "${2}" > "${1}"_"${2}"_mapped.sam;

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo "#######################################";
echo "#convert sam to bam#";
echo "#######################################";

samtools view -bS "${1}"_"${2}"_mapped.sam -o "${1}"_"${2}"_mapped.bam;

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo "#######################################";
echo "#indexing bam file#";
echo "#######################################";

samtools index "${1}"_"${2}"_mapped.bam;

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo "#######################################";
echo "# remove non-standard chromosomes#";
echo "#######################################";

sed '/chrM/d;/random/d;/chrUn/d' "${1}"_"${2}"_mapped.bam > "${1}"_"${2}"_mapped_remchr.bam
 
date +"%D %T":"##Process Completed Time##";
 
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo "#######################################";
echo "#sorting bam file#";
echo "#######################################";

samtools sort "${1}"_"${2}"_mapped_remchr.bam -o "${1}"_"${2}"_mapped_remchr_sorted.bam;

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo "#######################################";
echo "#remove duplicates#";
echo "#######################################";

samtools rmdup "${1}"_"${2}"_mapped_remchr_sorted.bam "${1}"_"${2}"_mapped_sorted_duprem.bam;

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo "#######################################";
echo "#add read groups to bam file#";
echo "#######################################";

java -jar "${gatk_root}"gatk.jar AddOrReplaceReadGroups I="${1}"_"${2}"_mapped_sorted_duprem.bam O="${1}"_"${2}"_mapped_sorted_duprem_RG.bam RGID=1 RGLB=lib2 RGPL=illumina RGPU=unit1 RGSM=3;

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo "#######################################";
echo "#indexing final bam file#";
echo "#######################################";

samtools index "${1}"_"${2}"_mapped_sorted_duprem_RG.bam;

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo "#######################################";
echo "#Base quality score Recalibration (BQSR)#";
echo "#######################################";

#Creating base recal table

java -jar "${gatk_root}"gatk.jar BaseRecalibrator -I "${1}"_"${2}"_mapped_sorted_duprem_RG.bam -R "${3}" --known-sites "${4}" -O base_recal.table;

date +"%D %T":"##Process Completed Time##";

#Apply BQSR
echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo "#######################################";
echo "#apply BQSR#";
echo "#######################################";

java -jar "${gatk_root}"gatk.jar ApplyBQSR -R "${3}" -I "${1}"_"${2}"_mapped_sorted_duprem_RG.bam --bqsr-recal-file base_recal.table -O "${1}"_"${2}"_mapped_sorted_duprem_RG_BQSR.bam;

date +"%D %T":"##Process Completed Time##";

echo -e "\n\n\n";
date +"%D %T":"##Process Start Time##";
echo "#######################################";
echo "#Variant calling using Haplotypecaller#";
echo "#######################################";

#haplotypecaller variant calling
java -jar "${gatk_root}"gatk.jar HaplotypeCaller -R "${3}" -I "${1}"_"${2}"_mapped_sorted_duprem_RG_BQSR.bam -O "${1}"_"${2}"_mapped_haplotypecaller.vcf;

date +"%D %T":"##Process Completed Time##";

echo done 

exit 0;
