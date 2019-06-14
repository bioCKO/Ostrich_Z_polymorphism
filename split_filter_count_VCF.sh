#!/usr/bin/bash
#______________________________________START CODE BLOCK 1_____________________________________#
# CODE BLOCK 1 deals with:
# Splitting vcf into Z scaffolds
# Filtering VCF
# Checking relatedness and genotypic sex
# Splitting vcf for male individuals
# Getting SNP counts
#_____________________________________________________________________________________________#
# Load the necessary modules
module load bioinfo-tools vcftools/0.1.15 bcftools/1.6 tabix plink/1.07
#*********************************************************************************************#
# VCF file
#*********************************************************************************************#
vcf='/proj/uppstore2017180/private/cornwallis/obsolete/data/interim/individual/all.bqsr.gatk.annotated.combined.snp.recal.vcf.gz'
#*********************************************************************************************#
# About the VCF:
#*********************************************************************************************#
# VCF contains data for three subspecies of ostriches: Black - Blue - Red
# black: P1878_107 - P1878_116, 
# blue: P1878_117 - P1878_126, 
# red: P1878_127 - P1878_136. 
# For each population, the first five individuals are male (e.g. P1878_107 -
# P1878_111 for black), the following five female. Note that one red female
# (P1878_133) has been misclassified and is rather a mix of black and
# blue. Moreover, P1878_122 is not clearly a blue sample. Both are removed from
# further analysis.
#*********************************************************************************************#
# Extract Z scaffolds from vcf
#*********************************************************************************************#
echo "Getting Z scaffolds"
vcftools --gzvcf ${vcf} --chr superscaffold26 --chr superscaffold54 --chr superscaffold35 --chr superscaffold36 --chr superscaffold62 --chr superscaffold67 --chr superscaffold69-1 --chr superscaffold93 --chr superscaffold63 --chr superscaffold88 --chr superscaffold83 --chr superscaffold92 --recode --recode-INFO-all --stdout > ../data/all.bqsr.gatk.annotated.snp.recal.Z.vcf
#*********************************************************************************************#
# bgzip and tabix the vcf file
bgzip ../data/all.bqsr.gatk.annotated.snp.recal.Z.vcf
tabix -p vcf ../data/all.bqsr.gatk.annotated.snp.recal.Z.vcf.gz
#*********************************************************************************************#
# Filtering criteria for SNPs
#*********************************************************************************************#
# 1. keep: keep the individuals of the subspecies given to this argument.
# 2. remove-filtered-all: Removes all sites wih a FILTER flag other than PASS.
# 3. remove-filtered-geno-all: Excludes all genotypes with a FILTER flag not equal to "." (a missing value) or PASS.
# 4. max-alleles
# 5. max-missing 1.0: No missing data is allowed.
#*********************************************************************************************#
echo "Split for subspecies and filter"
zvcf='all.bqsr.gatk.annotated.snp.recal.Z.vcf.gz'    
for subspecies in black.txt blue.txt red.txt # black.txt blue.txt red.txt are located in ostrich_Z_diversity/data/
do
species=$(echo $subspecies | cut -f1 -d ".")
vcfname=$(echo $zvcf | sed "s/all/$species/" | sed "s/.gz//")
echo $species $vcfname
echo "Filter VCF"
vcftools --gzvcf ../data/${zvcf} --keep ../data/${subspecies} --remove-filtered-all --remove-filtered-geno-all --max-alleles 2 --max-missing 1.0 --recode --recode-INFO-all --stdout > ../data/${vcfname}
bgzip ../data/${vcfname}
tabix -p vcf ../data/${vcfname}.gz
# Check relatedness among individuals
mkdir -p ../data/relatedness
vcftools --gzvcf ../data/${vcfname}.gz --relatedness2 --out ../data/relatedness/${subspecies}
# Check genotypic sex of the individuals
mkdir -p ../data/checkSex
vcftools --gzvcf ../data/${vcfname}.gz --chr superscaffold62 --chr superscaffold67 --chr superscaffold69-1 --chr superscaffold93 --chr superscaffold63 --chr superscaffold88 --chr superscaffold83 --chr superscaffold92 --recode --recode-INFO-all --stdout > ../data/checkSex/${species}.nonPAR.vcf
bgzip ../data/checkSex/${species}.nonPAR.vcf
tabix -p vcf ../data/checkSex/${species}.nonPAR.vcf.gz
vcftools --gzvcf ../data/checkSex/${species}.nonPAR.vcf.gz --plink --out ../data/checkSex/${species}.plink
awk '{print "X""\t"$2"\t"$3"\t"$4}' ../data/checkSex/${species}.plink.map > ../data/checkSex/${species}.plink.map2
rm -f ../data/checkSex/${species}.plink.map
mv ../data/checkSex/${species}.plink.map2 ../data/checkSex/${species}.plink.map 
plink --file ../data/checkSex/${species}.plink --check-sex --out ../data/checkSex/${species}.nonPAR
done
#*********************************************************************************************#
# hwe: filter all sites that deviate from HWE. This filters out sites that for example all
# male individuals are heterozygous.
echo "Split for males for each subspecies and filter"
zvcf='all.bqsr.gatk.annotated.snp.recal.Z.vcf.gz'    
for subspecies in black_male.txt blue_male.txt red_male.txt # black_male.txt blue_male.txt red_male.txt are located in ostrich_Z_diversity/data/
do
species=$(echo $subspecies | cut -f1 -d ".")
zvcf='all.bqsr.gatk.annotated.snp.recal.Z.vcf.gz'    
vcfname=$(echo $zvcf | sed "s/all/$species/" | sed "s/.gz//")
echo $species $vcfname
echo "Filter VCF"
# Here I add "maf" filtering of 0.1, because the subspecies specific vcf may contain sites where
# it is variable in 2 subspecies and not in 1. In that case, the species will have a site 
# with frequency of REF 10 and ALT 0. The site may also be a substituition compared to the 
# reference. These cases should be less frequent and may mostly represent errors.
vcftools --gzvcf ../data/${zvcf} --keep ../data/${subspecies} --remove-filtered-all --remove-filtered-geno-all --max-alleles 2 --max-missing 1.0 --maf 0.05 --hwe 0.05 --recode --recode-INFO-all --stdout > ../data/${vcfname}
bgzip ../data/${vcfname}
tabix -p vcf ../data/${vcfname}.gz
done
#*********************************************************************************************#
echo "Take allele counts of each subspecies"
for species in black blue red
do
vcftools --gzvcf ../data/${species}_male.bqsr.gatk.annotated.snp.recal.Z.vcf.gz --counts --out ../data/${species}.male.counts
done
#________________________________________END CODE BLOCK 1_____________________________________#
#___________________________________________SCRIPT NAME_______________________________________#
#___________________________________split_filter_count_VCF.sh_________________________________#
#_____________________________________________RUN_____________________________________________#
#_________________________________bash split_filter_count_VCF.sh______________________________#
#_____________________________________________________________________________________________#
#_____________________________________________________________________________________________#
#_____________________________________________________________________________________________#
