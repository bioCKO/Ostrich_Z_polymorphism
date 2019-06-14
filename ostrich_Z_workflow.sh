#*********************************************************************************************#
#*********************************************************************************************#
#****************************** Analysis Steps for the project *******************************#
#************************ Genetic diversity on ostrich Z chromosome **************************#
#************************************* Homa Papoli Yazdi *************************************#
#*********************************************************************************************#
#*********************************************************************************************#

#*********************************************************************************************#
# The working directory is located in:
#*********************************************************************************************#
# /proj/uppstore2017180/private/homap/ostrich_Z_diversity
#*********************************************************************************************#
#************************************** ATTENTION ********************************************#
#************************ UPPSTORE PROJECTS ARE NOT BACKED UP ********************************# 
#*********************************************************************************************#
# The final directory is also copied in /home/homap on rackham.
# The final directory is also copied into Dropbox account: hpapoli@gmail.com.
#*********************************************************************************************#
# Structure of the working directory located in /proj/uppstore2017180/private/homap
# is as follows:
# ostrich_Z_diversity/
# ├── data
# ├── result
# │   ├── graphs
# │   └── tables
# └── scripts
# NB: Structure was produced by command 'tree' at the start of the project.
#*********************************************************************************************#
#**************************** ALL SCRIPTS ARE RUN FROM /scripts ******************************#
#*********************************************************************************************#    
# In this document, the code in all bash scripts used is included. At the end of each block of 
# code, the script name containing that code block and how it was run is written. In practice,
# this script can be used to reproduce all results of the project. The necessary python scripts
# are found in the /scripts directory. 
#*********************************************************************************************#
#*********************************************************************************************#
#*********************************************************************************************#
#*********************************************************************************************#
# PART I: 
#		A. Calculating diversity (pi) in black, blue and red subspecies across Z
#		B. Calculate the density of GC, gene, repeat and get recombination rate, genetic 
#		   distance and physical distance from the pseudo-autosomal boundary
#		C. Produce graphs of diversity in the 3 subspecies across Z
#		D. Produce graphs of diversity in the 3 subspecies across PAR boundary
#		E. Use diversity (pi) as the response variable and the factors in (B) as the 
#		   explanatory variables and conduct multiple linear regression using GLM
#*********************************************************************************************#
#*********************************************************************************************#
#*********************************************************************************************#
#*********************************************************************************************#
# The response variable is Z chromosome heterozygosity: 
#
# Response _______ SNP diversity: Neutral genetic diversity measured by intergenic sequences 
#                                 at least 1000 bp away from genes.
#*********************************************************************************************#
# The explanatory variables for levels of heterozygosity are the following:
#
#             _____ Recombination rate: The genetic map obtained in Yazdi & Ellegren 2018 is used.
#            |                          Sex average, male and female maps were used to calculate                                          
#            |                          Recombination rate in cM/Mb.
#            |                          
# Explanatory _____ GC content 
# 		     |
# 		     |
# 		     |_____ Gene density
#            |
#            |
#            |_____ Repeat density
#            |
#            |
#            |_____ Genetic distance from the pseudo-autosomal (PAR) boundary
#            |
#            |
#            |_____ Physical distance from the PAR boundary
#
#*********************************************************************************************#
# Genetic diversity analysis strategy:
# 
# Window based analysis: Window sizes 100 Kb - 200 Kb - 500 Kb - 1 Mb
# Each value of response and explanatory is measured within one window and a mean over the
# window is taken. 
# Genetic diversity (pi) is calculated using the number of SNPs within the intergenic regions 
# of each window at least 1000 bp away from the genes. To get pi per site, average pi per
# window is divided by the total number of intergenic bases that have passed the filtering
# criteria: minimum coverage + masked for repeats
#
#                        	     Window 1
#
# Ind1 |ATG[C/T]GTGTGGNNNNNNNNATGGTGC[A/G]CGTGACGTCTNNNNNNNNNNNNNNNNTGTACACGTG[G/C]GTGACACACGT| 
# Ind2 |ATG[C/C]GTGTGGNNNNNNNNATGGTGC[A/G]CGTGACGTCTNNNNNNNNNNNNNNNNTGTACACGTG[G/C]GTGACACACGT|
# Ind3 |ATG[C/C]GTGTGGNNNNNNNNATGGTGC[A/A]CGTGACGTCTNNNNNNNNNNNNNNNNTGTACACGTG[G/G]GTGACACACGT|
# Ind4 |ATG[C/C]GTGTGGNNNNNNNNATGGTGC[A/G]CGTGACGTCTNNNNNNNNNNNNNNNNTGTACACGTG[C/C]GTGACACACGT|
# Ind5 |ATG[C/C]GTGTGGNNNNNNNNATGGTGC[A/A]CGTGACGTCTNNNNNNNNNNNNNNNNTGTACACGTG[C/C]GTGACACACGT|
#                     Repeat                        Coverage filter                        
#																	
# To obtain confidence interval of the mean of pi per window, we perform resampling of all 
# intergenic regions per window for 1000 times.
#*********************************************************************************************#
#*********************************************************************************************#
#**************************************** CODE BLOCKS ****************************************#
#*********************************************************************************************#
#*********************************************************************************************#
#_____________________________________________________________________________________________#
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


#!/usr/bin/bash
#______________________________________START CODE BLOCK 2_____________________________________#
# Load the necessary modules
module load bioinfo-tools BEDTools
#*********************************************************************************************#
# CODE BLOCK 2 deals with:
# Correcting GFF of ostrich genome assembly
# Obtaining coordinates of intergenic, intronic and CDS regions of GFF in BED format
#*************************************************************************************************************************************************************************************************************************
# CDS
#*************************************************************************************************************************************************************************************************************************
# Correct the gff format
#*************************************************************************************************************************************************************************************************************************
echo "Reformatting GFF"
./reformat_gff.py > ../data/Struthio_camelus.OM.gene.20130116.reformated.gff
#*************************************************************************************************************************************************************************************************************************
# Remove overlapping mRNA
#*************************************************************************************************************************************************************************************************************************
echo "Removing overlapping mRNA"
awk '$3=="mRNA"' ../data/Struthio_camelus.OM.gene.20130116.reformated.gff | awk '{print $1"\t"$4-1"\t"$5}' | sort -k1,1 -k2,2n | bedtools merge -i - -c 1 -o count | awk '{if($4==1) print $0}' | bedtools intersect -a ../data/Struthio_camelus.OM.gene.20130116.reformated.gff -b - > ../data/Struthio_camelus.OM.gene.20130116.reformated.uniq.gff
#*************************************************************************************************************************************************************************************************************************
# Extract CDS coordinates in bed format
#*************************************************************************************************************************************************************************************************************************
echo "Extract CDS coordinates"
./extract_from_gff.py ../data/Struthio_camelus.OM.gene.20130116.reformated.uniq.gff CDS all > ../data/ostrich.all.CDS.coord
#*************************************************************************************************************************************************************************************************************************
# Extract intronic coordinates
#*************************************************************************************************************************************************************************************************************************
echo "Extract intronic coordinates"
./extract_from_gff.py ../data/Struthio_camelus.OM.gene.20130116.reformated.uniq.gff intron all > ../data/ostrich.all.intron.coord
#*************************************************************************************************************************************************************************************************************************
# Extract intergenic coordinates
#*************************************************************************************************************************************************************************************************************************
echo "Extract intergenic coordinates"
sort -k1,1 -k2,2n ../data/Struthio_camelus.20130116.OM.length.txt > ../data/Struthio_camelus.20130116.OM.length.sorted.txt
awk '$3=="mRNA"' ../data/Struthio_camelus.OM.gene.20130116.reformated.uniq.gff | awk '{print $1"\t"$4-1"\t"$5}' | sort -k1,1 -k2,2n | bedtools merge -i - -c 1 -o count | bedtools complement -i - -g ../data/Struthio_camelus.20130116.OM.length.sorted.txt > ../data/ostrich.intergene.coord
#________________________________________END CODE BLOCK 2_____________________________________#
#___________________________________________SCRIPT NAME_______________________________________#
#________________________________extract_ostrich_coordinates.sh_______________________________#
#_____________________________________________RUN_____________________________________________#
#______________________________bash extract_ostrich_coordinates.sh____________________________#
#_____________________________________________________________________________________________#
#_____________________________________________________________________________________________#


#!/usr/bin/bash
#______________________________________START CODE BLOCK 3_____________________________________#
# CODE BLOCK 2 deals with:
# Background DNA must be masked for repeats and coverage.
# |ATG[C/C]GTGTGGNNNNNNNNATGGTGC[A/A]CGTGACGTCTNNNNNNNNNNNNNNNNTGTACACGTG[C/C]GTGACACAACGT|
#                 Repeat                        Coverage filter
# Hard mask background for repeats and coverage. Use subspecies-specific background coverage.
# Repeats are already hard masked using the repeat masker. Look in the RepeatMasker_Ostrich.sh  
#_____________________________________________________________________________________________#
# Load the necessary modules
module load bioinfo-tools BEDTools
#*********************************************************************************************#
# Background coverage
#*********************************************************************************************#
echo "Hard masking fasta for coverage"

coverageAddr='/proj/uppstore2017180/private/cornwallis/data/interim/individual/genomecov'

for subspecies in black blue red
do
echo ${subspecies}
echo "Sorting the coverage bed file"
zcat ${coverageAddr}/${subspecies}.coverage.7X.70pct.bed.gz | sort -k1,1 -k2n,2 > ../data/${subspecies}.coverage.7X.70pct.sorted.bed
echo "Return complement bed file with bases having low coverage"
bedtools complement -i ../data/${subspecies}.coverage.7X.70pct.sorted.bed -g ../data/Struthio_camelus.20130116.OM.length.sorted.txt > ../data/${subspecies}.coverage.7X.70pct.sorted.complement.bed
echo "Masking fasta"
bedtools maskfasta -fi ../data/Ostrich_repeatMask/Struthio_camelus.20130116.OM.fa.masked -bed ../data/${subspecies}.coverage.7X.70pct.sorted.complement.bed -fo ../data/${subspecies}.repeat.cov.masked.fa
done
#________________________________________END CODE BLOCK 3_____________________________________#
#___________________________________________SCRIPT NAME_______________________________________#
#____________________________________hard_masking_background.sh_______________________________#
#_____________________________________________RUN_____________________________________________#
#_________________________________bash hard_masking_background.sh_____________________________#
#_____________________________________________________________________________________________#
#_____________________________________________________________________________________________#


#!/usr/bin/bash
#______________________________________START CODE BLOCK 4_____________________________________#
# CODE BLOCK 4 deals with:
# Creating windows of 100, 200, 500 and 1000 Kb
# Get an overlap of each window file with repeats, CDS, intronic and intergenic sequences
# Do the correction if the same feature coordinate falls into two windows. That is:
# 100 200 50 210 160
# 200 300 50 210 160
# Get the sum of overlapping features per window
#_____________________________________________________________________________________________
# Load the necessary modules
module load bioinfo-tools BEDTools
#*********************************************************************************************#
# Get windows per scaffolds
#*********************************************************************************************#
echo "Create windows"
for window in 100Kb 200Kb 500Kb 1000Kb
do
w_size=$(echo $window | sed 's/Kb/000/g')
python sliding_window_ZA.py ../data/Z_scaffolds_length.txt ${w_size} > ../data/ostrich.Z.${window}.bed
done
#*********************************************************************************************#
# Get overlap between windows and element
#*********************************************************************************************#
# Repeats
echo "Repeats"
#*********************************************************************************************#
# merge overlapping repeat coordinates
repeat_addr='../data/Ostrich_repeatMask/Struthio_camelus.20130116.OM.fa.out'
awk -v OFS="\t" '$1=$1' ${repeat_addr} | awk '{$6=$6-1; print $5"\t"$6"\t"$7}' | awk 'NR>2' | mergeBed -i - | grep -f ../data/Z_scaffolds.txt > ../data/SC.allRepeats.bed
for window in 100Kb 200Kb 500Kb 1000Kb
do
echo ${window}
bedtools intersect -a ../data/ostrich.Z.${window}.bed -b ../data/SC.allRepeats.bed -wao > ../data/ostrich.Z.${window}.AllRepeat.overlap.txt
python get_sum_overlapping.py ../data/ostrich.Z.${window}.AllRepeat.overlap.txt > ../data/ostrich.Z.${window}.AllRepeat.overlap.sum.txt
# Use an already repeat masked genome
# Because the repeat masked genome is used, sum of bases and gc for repeat elements is zero in the outout
# as due to hard masking, these positions are set to "N" in the repeat masked fasta file. 
for subspecies in black blue red
do
echo ${subspecies}
fasta=../data/${subspecies}.repeat.cov.masked.fa
python get_density_feature.py ../data/ostrich.Z.${window}.AllRepeat.overlap.sum.txt ${fasta} ../data/${subspecies}.Z.${window}.AllRepeat.overlap.density.txt 
sort -k1,1 -k2n,2 ../data/${subspecies}.Z.${window}.AllRepeat.overlap.density.txt | awk 'NR>1' > ../data/${subspecies}.Z.${window}.AllRepeat.overlap.density.sorted.txt
done
done
#*************************************************************************************************************************************************************************************************************************
# CDS
echo "CDS"
#*************************************************************************************************************************************************************************************************************************
for window in 100Kb 200Kb 500Kb 1000Kb
do
echo $window
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/ostrich.all.CDS.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | bedtools intersect -a ../data/ostrich.Z.${window}.bed -b - -wao > ../data/ostrich.Z.${window}.CDS.overlap.txt	
python get_sum_overlapping.py ../data/ostrich.Z.${window}.CDS.overlap.txt	 > ../data/ostrich.Z.${window}.CDS.overlap.sum.txt
for subspecies in black blue red
do
echo ${subspecies}
fasta=../data/${subspecies}.repeat.cov.masked.fa
python get_density_feature.py ../data/ostrich.Z.${window}.CDS.overlap.sum.txt ${fasta} ../data/${subspecies}.Z.${window}.CDS.overlap.density.txt
sort -k1,1 -k2n,2 ../data/${subspecies}.Z.${window}.CDS.overlap.density.txt | awk 'NR>1' > ../data/${subspecies}.Z.${window}.CDS.overlap.density.sorted.txt
done
done
#*************************************************************************************************************************************************************************************************************************
# Intron
echo "Intron"
#*************************************************************************************************************************************************************************************************************************
for window in 100Kb 200Kb 500Kb 1000Kb
do
echo ${window}
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/ostrich.all.intron.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | bedtools intersect -a ../data/ostrich.Z.${window}.bed -b - -wao > ../data/ostrich.Z.${window}.intronic.overlap
python get_sum_overlapping.py ../data/ostrich.Z.${window}.intronic.overlap > ../data/ostrich.Z.${window}.intronic.overlap.sum.txt
for subspecies in black blue red
do
echo ${subspecies}	
fasta=../data/${subspecies}.repeat.cov.masked.fa
python get_density_feature.py ../data/ostrich.Z.${window}.intronic.overlap.sum.txt ${fasta} ../data/${subspecies}.Z.${window}.intronic.overlap.density.txt
sort -k1,1 -k2n,2 ../data/${subspecies}.Z.${window}.intronic.overlap.density.txt | awk 'NR>1' > ../data/${subspecies}.Z.${window}.intronic.overlap.density.sorted.txt
done
done
#*************************************************************************************************************************************************************************************************************************
# Intergenic
echo "Intergenic"
#*************************************************************************************************************************************************************************************************************************
for window in 100Kb 200Kb 500Kb 1000Kb
do
bedtools intersect -a ../data/ostrich.Z.${window}.bed -b ../data/ostrich.intergene.coord -wao > ../data/ostrich.Z.${window}.intergenic.overlap
python get_sum_overlapping.py ../data/ostrich.Z.${window}.intergenic.overlap > ../data/ostrich.Z.${window}.intergenic.overlap.sum.txt
for subspecies in black blue red
do
echo ${subspecies}	
fasta=../data/${subspecies}.repeat.cov.masked.fa
python get_density_feature.py ../data/ostrich.Z.${window}.intergenic.overlap.sum.txt ${fasta} ../data/${subspecies}.Z.${window}.intergenic.overlap.density.txt
sort -k1,1 -k2n,2 ../data/${subspecies}.Z.${window}.intergenic.overlap.density.txt | awk 'NR>1' > ../data/${subspecies}.Z.${window}.intergenic.overlap.density.sorted.txt
done
done
#________________________________________END CODE BLOCK 4_____________________________________#
#___________________________________________SCRIPT NAME_______________________________________#
#____________________________________window_feature_overlap.sh________________________________#
#_____________________________________________RUN_____________________________________________#
#________________________________bash window_feature_overlap.sh_______________________________#
#_____________________________________________________________________________________________#
#_____________________________________________________________________________________________#


#!/usr/bin/bash
#______________________________________START CODE BLOCK 5_____________________________________#
# CODE BLOCK 5 deals with:
# SNP annotations for low coverage and repeats (as done for background) plus annotation for
# intergenic, intronic and CDS categories.
#_____________________________________________________________________________________________#
# Load the necessary modules
module load bioinfo-tools BEDTools
#*********************************************************************************************#
# Remove the first 1Kb and the last 1Kb of intergenic coordinates. This achieves 3 goals:
# 1. Removes the first and last 1Kb of scaffolds, prone to errors. 
# 2. Removes the 1Kb closest to the genes, hence, reducing the effect of linked selection.
# 3. Removes very short scaffolds or contigs of =< 1000 bps.
awk '{print $1"\t"$2+1000"\t"$3-1000}' ../data/ostrich.intergene.coord | awk '$3>0' | awk '$3>$2' > ../data/ostrich.intergene.1Kb.removed.coord
for subspecies in black blue red
do
echo ${subspecies} 
# Annotate each SNP for low coverage
awk '{if(NR>1) print $1"\t"$2-1"\t"$2"\t"$4"\t"$5"\t"$6}' ../data/${subspecies}.male.counts.frq.count | bedtools intersect -a - -b ../data/${subspecies}.coverage.7X.70pct.sorted.complement.bed -wao | awk '{if($7==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""PASS"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""LOW_COV"}' > ../data/${subspecies}.male.counts.lowcov
# Annotate each SNP for repeat 
bedtools intersect -a ../data/${subspecies}.male.counts.lowcov -b ../data/SC.allRepeats.bed -wao | awk '{if($8==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""PASS"; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""REPEAT"}' > ../data/${subspecies}.male.counts.repeat
# Annotate each SNP for intergenic
bedtools intersect -a ../data/${subspecies}.male.counts.repeat -b ../data/ostrich.intergene.1Kb.removed.coord -wao | awk '{if($9==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t""."; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t""INTERGENIC"}' | sort -k1,1 -k2n,2 > ../data/${subspecies}.male.counts.repeat.intergenic
# Annotate each SNP for intronic
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/ostrich.all.intron.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | bedtools intersect -a ../data/${subspecies}.male.counts.repeat.intergenic -b - -wao | awk '{if($10==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t""."; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t""INTRON"}' > ../data/${subspecies}.male.counts.repeat.intergenic.intronic
# Annotate each SNP for CDS 
awk '{if($2>$3) print $1"\t"$3"\t"$2; else print $0}' ../data/ostrich.all.CDS.coord | awk '{split($1, a, "_"); print a[1]"\t"$2"\t"$3}' | bedtools intersect -a ../data/${subspecies}.male.counts.repeat.intergenic.intronic -b - -wao | awk '{if($11==".") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t""."; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t""CDS"}' > ../data/${subspecies}.male.counts.repeat.intergenic.intronic.cds
done
#*************************************************************************************************************************************************************************************************************************
for window in 100Kb 200Kb 500Kb 1000Kb
do
echo ${window}
for subspecies in black blue red
do
echo ${subspecies}
bedtools intersect -a ../data/${subspecies}.male.counts.repeat.intergenic.intronic.cds -b ../data/${subspecies}.Z.${window}.AllRepeat.overlap.density.sorted.txt -wao | awk '{$12=$21=""; print $0}' | sed 's/  / /g' | sed 's/ /\t/g' | sed 's/\t$//' > ../data/${subspecies}.${window}.repeat.txt
bedtools intersect -a ../data/${subspecies}.${window}.repeat.txt -b ../data/${subspecies}.Z.${window}.CDS.overlap.density.sorted.txt -wao | awk '{$20=$21=$22=$23=$24=$29=""; print $0}' | sed 's/      / /g' | sed 's/  / /g' | sed 's/ /\t/g' | sed 's/\t$//' > ../data/${subspecies}.${window}.repeat.CDS.txt
bedtools intersect -a ../data/${subspecies}.${window}.repeat.CDS.txt -b ../data/${subspecies}.Z.${window}.intergenic.overlap.density.sorted.txt -wao | awk '{$24=$25=$26=$27=$28=$33=""; print $0}' | sed 's/      / /g' | sed 's/  / /g' | sed 's/ /\t/g' | sed 's/\t$//' > ../data/${subspecies}.${window}.repeat.CDS.intergenic.txt
bedtools intersect -a ../data/${subspecies}.${window}.repeat.CDS.intergenic.txt -b ../data/${subspecies}.Z.${window}.intronic.overlap.density.sorted.txt -wao | awk '{$28=$29=$30=$31=$32=$37=""; print $0}' | sed 's/      / /g' | sed 's/  / /g' | sed 's/ /\t/g' | sed 's/\t$//' > ../data/${subspecies}.${window}.repeat.CDS.intergenic.intronic.txt
cat ../data/header ../data/${subspecies}.${window}.repeat.CDS.intergenic.intronic.txt > ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.txt
done
done
#________________________________________END CODE BLOCK 5_____________________________________#
#___________________________________________SCRIPT NAME_______________________________________#
#________________________________________SNP_annotation.sh____________________________________#
#_____________________________________________RUN_____________________________________________#
#______________________________________bash SNP_annotation.sh_________________________________#
#_____________________________________________________________________________________________#
#_____________________________________________________________________________________________#


#!/usr/bin/bash
#______________________________________START CODE BLOCK 6_____________________________________#
# Calculate pi per site
for window in 100Kb 200Kb 500Kb 1000Kb
do
echo ${window}
for subspecies in black blue red
do
./pi_dxy.py ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.txt ../data/${subspecies}.${window}.success.failure.txt #> ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.txt
grep -v "REPEAT" ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.txt | grep -v "LOW_COV" > ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.norepeat.txt
grep -v "REPEAT" ../data/${subspecies}.${window}.success.failure.txt | grep -v "LOW_COV" > ../data/${subspecies}.${window}.success.failure.norepeat.txt
done
done
#*************************************************************************************************************************************************************************************************************************
# Separate pi into functional categories and remove repeat elements
for window in 100Kb 200Kb 500Kb 1000Kb
do
echo ${window}
for subspecies in black blue red
do
echo ${subspecies}
for element in INTERGENIC INTRON CDS
do
echo ${element}
grep ${element} ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.norepeat.txt > ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.${element}.txt
grep ${element} ../data/${subspecies}.${window}.success.failure.norepeat.txt > ../data/${subspecies}.${window}.success.failure.norepeat.${element}.txt
done
done
done
#*************************************************************************************************************************************************************************************************************************
# Get sum of pi per window
for window in 100Kb 200Kb 500Kb 1000Kb
do
echo ${window}
for subspecies in black blue red
do
echo ${subspecies}
python sum_pi_window_intergenic.py ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.INTERGENIC.txt > ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.INTERGENIC.sum.txt
python sum_pi_window_intronic.py ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.INTRON.txt > ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.INTRON.sum.txt
python sum_pi_window_CDS.py ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.CDS.txt > ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.CDS.sum.txt
python sum_pi_window.py ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.norepeat.txt > ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.norepeat.sum.txt
python sum_pi_window_intergenic_2.py ../data/${subspecies}.${window}.success.failure.norepeat.INTERGENIC.txt > ../data/${subspecies}.${window}.success.failure.norepeat.INTERGENIC.sum.txt 
done
done
#________________________________________END CODE BLOCK 6_____________________________________#
#___________________________________________SCRIPT NAME_______________________________________#
#_________________________________________calculate_pi.sh_____________________________________#
#_____________________________________________RUN_____________________________________________#
#_______________________________________bash calculate_pi.sh__________________________________#
#_____________________________________________________________________________________________#
#_____________________________________________________________________________________________#


#!/usr/bin/bash
#______________________________________START CODE BLOCK 7_____________________________________#
# Add confidence interval for intergenic regions
for window in 200Kb 500Kb 1000Kb
do
echo ${window}
for subspecies in black blue red
do
echo ${subspecies}
python resampling_scaffold.py ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.INTERGENIC.txt > ../data/${subspecies}.${window}.INTERGENIC.mean.CI.txt
python scaffold_to_chr_vcf.py ../data/${subspecies}.${window}.INTERGENIC.mean.CI.txt > ../data/${subspecies}.${window}.INTERGENIC.mean.CI.ChrZ.txt
done
done
#________________________________________END CODE BLOCK 7_____________________________________#
#___________________________________________SCRIPT NAME_______________________________________#
#____________________________________resampling_intergenic.sh_________________________________#
#_____________________________________________RUN_____________________________________________#
#_________________________________bash resampling_intergenic.sh_______________________________#
#_____________________________________________________________________________________________#
#_____________________________________________________________________________________________#


#!/usr/bin/bash
#______________________________________START CODE BLOCK 8_____________________________________#
# Obtain per window sex averaged, male and female recombination rates
#_____________________________________________________________________________________________#
# Load the necessary modules
module load bioinfo-tools BEDTools
#_____________________________________________________________________________________________#
./recombination.py ../data/LGZ3.sex_averaged.lifted.bed > ../data/LGZ3.sex_averaged.lifted.rec.rate.txt
sed 's/superscaffold54.1/superscaffold54/g' ../data/LGZ3.sex_averaged.lifted.rec.rate.txt > ../data/LGZ3.sex_averaged.lifted.rec.rate.54.txt
awk '{if($2>$3) print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/LGZ3.sex_averaged.lifted.rec.rate.54.txt > ../data/LGZ3.sex_averaged.lifted.rec.rate.54.corrected.txt

for window in 200Kb 500Kb 1000Kb
do
echo ${window}
for subspecies in black blue red
do 
echo ${subspecies}
bedtools intersect -a ../data/ostrich.Z.${window}.bed -b ../data/LGZ3.sex_averaged.lifted.rec.rate.54.corrected.txt -wao > ../data/LGZ3.sex_averaged.lifted.rec.rate.${window}.txt
./get_recombination_per_window.py ../data/LGZ3.sex_averaged.lifted.rec.rate.${window}.txt ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.INTERGENIC.sum.txt > ../data/${subspecies}_rec_per_${window}.txt
done
done

for window in 200Kb 500Kb 1000Kb
do
echo ${window}
for subspecies in black blue red
do 
echo ${subspecies}
cut -f4 ../data/${subspecies}_rec_per_${window}.txt | paste ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.INTERGENIC.sum.txt - > ../data/${subspecies}.${window}.INTERGENIC.rec.sum.txt 
awk 'BEGIN{print "ChrZ_start""\t""ChrZ_end""\t""Resampling_mean""\t""Low_CI""\t""Up_CI""\t""SD"}{print $5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' ../data/${subspecies}.${window}.INTERGENIC.mean.CI.ChrZ.txt | paste ../data/${subspecies}.${window}.INTERGENIC.rec.sum.txt - > ../data/${subspecies}.${window}.INTERGENIC.rec.resampling.Z.sum.txt
done
done
#*********************************************************************************************#
# Get male and female recombination rate for use in multiple linear regression model and SLIM simulation
grep 'female_input' ../data/LGZ3.male_female.lifted.bed > ../data/LGZ3.female.tmp
grep -w 'male_input' ../data/LGZ3.male_female.lifted.bed > ../data/LGZ3.male.tmp
./recombination.py ../data/LGZ3.female.tmp > ../data/LGZ3.female.tmp.rec.rate.txt
./recombination.py ../data/LGZ3.male.tmp > ../data/LGZ3.male.tmp.rec.rate.txt

sed 's/superscaffold54.1/superscaffold54/g' ../data/LGZ3.female.tmp.rec.rate.txt > ../data/LGZ3.female.tmp.rec.rate.54.txt
sed 's/superscaffold54.1/superscaffold54/g' ../data/LGZ3.male.tmp.rec.rate.txt > ../data/LGZ3.male.tmp.rec.rate.54.txt

awk '{if($2>$3) print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/LGZ3.female.tmp.rec.rate.54.txt > ../data/LGZ3.female.tmp.rec.rate.54.corrected.txt
awk '{if($2>$3) print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ../data/LGZ3.male.tmp.rec.rate.54.txt > ../data/LGZ3.male.tmp.rec.rate.54.corrected.txt

for window in 200Kb 500Kb 1000Kb
do
echo ${window}
for sex in male female
do
for subspecies in black blue red
do
echo ${subspecies} ${sex}
bedtools intersect -a ../data/ostrich.Z.${window}.bed -b ../data/LGZ3.${sex}.tmp.rec.rate.54.corrected.txt -wao > ../data/LGZ3.${sex}.lifted.rec.rate.window_overlap.${window}.txt
./get_recombination_per_window.py ../data/LGZ3.${sex}.lifted.rec.rate.window_overlap.${window}.txt ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.INTERGENIC.sum.txt > ../data/${subspecies}_rec_${sex}_per_${window}.txt
done
done
done

for window in 200Kb 500Kb 1000Kb
do
echo ${window}
for subspecies in black blue red
do
echo ${subspecies}	
cut -f4 ../data/${subspecies}_rec_male_per_${window}.txt | awk 'BEGIN{print "male_CM_per_bp"}{if(NR>1) print $1}' | paste ../data/${subspecies}.${window}.INTERGENIC.rec.resampling.Z.sum.txt - > ../data/${subspecies}.${window}.INTERGENIC.rec.resampling.Z.sex.txt   
cut -f4 ../data/${subspecies}_rec_female_per_${window}.txt | awk 'BEGIN{print "female_CM_per_bp"}{if(NR>1) print $1}' | paste ../data/${subspecies}.${window}.INTERGENIC.rec.resampling.Z.sex.txt - >  ../data/${subspecies}.${window}.INTERGENIC.rec.resampling.Z.sex2.txt   
done
done
rm -f ../data/${subspecies}.${window}.INTERGENIC.rec.resampling.Z.sex.txt 
#*************************************************************************************************************************************************************************************************************************

#________________________________________END CODE BLOCK 8_____________________________________#
#___________________________________________SCRIPT NAME_______________________________________#
#_______________________________________recombinaton_rate.sh__________________________________#
#_____________________________________________RUN_____________________________________________#
#___________________________________bash recombinaton_rate.sh_________________________________#
#_____________________________________________________________________________________________#
#_____________________________________________________________________________________________#

#!/usr/bin/bash
#______________________________________START CODE BLOCK 9_____________________________________#
# Obtain genetic distance from the centre of a window to the PAR boundary
#_____________________________________________________________________________________________#
# Load the necessary modules
module load bioinfo-tools BEDTools
#_____________________________________________________________________________________________#
for window in 200Kb 500Kb 1000Kb
do
echo ${window}
awk '{print "ChrZ""\t"$9"\t"$10"\t"$8}' ../data/black.${window}.INTERGENIC.rec.resampling.Z.sex2.txt | awk 'NR>1' | awk '{if($2>$3) print $1"\t"$3"\t"$2"\t"$4; else print $1"\t"$2"\t"$3"\t"$4}' > ../data/black.${window}.INTERGENIC.rec.resampling.Z.sex3.txt
python CM_to_PAR_boundary.py  ../data/black.${window}.INTERGENIC.rec.resampling.Z.sex2.txt > ../data/black.${window}.INTERGENIC.rec.resampling.Z.sex4.txt
bedtools intersect -a ../data/black.${window}.INTERGENIC.rec.resampling.Z.sex4.txt -b ../data/black.${window}.INTERGENIC.rec.resampling.Z.sex3.txt -wao > ../data/black.${window}.INTERGENIC.rec.overlap.txt
awk '{print $1"\t"$5"\t"$6}' ../data/black.${window}.INTERGENIC.mean.CI.ChrZ.txt > test4
python CM_to_PAR_boundary_step2.py ../data/black.${window}.INTERGENIC.rec.overlap.txt test4 > ../data/distance_from_PAR_${window}.txt
awk 'NR>1' ../data/distance_from_PAR_${window}.txt | sort -k2n,2 | awk 'BEGIN{print "Scaffold""\t""Focal_site""\t""PAR_boundary""\t""CM_from_PAR_boundary"}{print $0}' > ../data/distance_from_PAR_${window}.sorted.txt
cut -f2,3,4 ../data/distance_from_PAR_${window}.sorted.txt | paste ../data/black.${window}.INTERGENIC.rec.resampling.Z.sex2.txt - > ../data/black.${window}.INTERGENIC.rec.resampling.Z.sex2.CMPAR.txt
cut -f2,3,4 ../data/distance_from_PAR_${window}.sorted.txt | paste ../data/blue.${window}.INTERGENIC.rec.resampling.Z.sex2.txt - > ../data/blue.${window}.INTERGENIC.rec.resampling.Z.sex2.CMPAR.txt
cut -f2,3,4 ../data/distance_from_PAR_${window}.sorted.txt | paste ../data/red.${window}.INTERGENIC.rec.resampling.Z.sex2.txt - > ../data/red.${window}.INTERGENIC.rec.resampling.Z.sex2.CMPAR.txt
done
#________________________________________END CODE BLOCK 9_____________________________________#
#___________________________________________SCRIPT NAME_______________________________________#
#____________________________________CM_from_PAR_boundary.sh__________________________________#
#_____________________________________________RUN_____________________________________________#
#__________________________________bash CM_from_PAR_boundary.sh_______________________________#
#_____________________________________________________________________________________________#
#_____________________________________________________________________________________________#

#!/usr/bin/bash
#______________________________________START CODE BLOCK 10____________________________________#
# Print the columns with an easier-to-read order for statistical analysis
for window in 200Kb 500Kb 1000Kb
do
echo ${window}
awk '{print $4"\t"$5}' ../data/black.${window}.success.failure.norepeat.INTERGENIC.sum.txt > ../data/black.${window}.succ.fail
awk '{print $4"\t"$5}' ../data/blue.${window}.success.failure.norepeat.INTERGENIC.sum.txt > ../data/blue.${window}.succ.fail
awk '{print $4"\t"$5}' ../data/red.${window}.success.failure.norepeat.INTERGENIC.sum.txt > ../data/red.${window}.succ.fail
awk '{if(NR>1 && $9>$10) print $1"\t"$2"\t"$3"\t"$10"\t"$9"\t"$4"\t"$11"\t"$12"\t"$13"\t"$14"\t"$8"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$5"\t"$6"\t"$7; else print $1"\t"$2"\t"$3"\t"$9"\t"$10"\t"$4"\t"$11"\t"$12"\t"$13"\t"$14"\t"$8"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$5"\t"$6"\t"$7}' ../data/black.${window}.INTERGENIC.rec.resampling.Z.sex2.CMPAR.txt | paste - ../data/black.${window}.succ.fail > ../data/black.${window}.INTERGENIC.rec.resampling.Z.forR.txt
awk '{if(NR>1 && $9>$10) print $1"\t"$2"\t"$3"\t"$10"\t"$9"\t"$4"\t"$11"\t"$12"\t"$13"\t"$14"\t"$8"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$5"\t"$6"\t"$7; else print $1"\t"$2"\t"$3"\t"$9"\t"$10"\t"$4"\t"$11"\t"$12"\t"$13"\t"$14"\t"$8"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$5"\t"$6"\t"$7}' ../data/blue.${window}.INTERGENIC.rec.resampling.Z.sex2.CMPAR.txt | paste - ../data/blue.${window}.succ.fail > ../data/blue.${window}.INTERGENIC.rec.resampling.Z.forR.txt
awk '{if(NR>1 && $9>$10) print $1"\t"$2"\t"$3"\t"$10"\t"$9"\t"$4"\t"$11"\t"$12"\t"$13"\t"$14"\t"$8"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$5"\t"$6"\t"$7; else print $1"\t"$2"\t"$3"\t"$9"\t"$10"\t"$4"\t"$11"\t"$12"\t"$13"\t"$14"\t"$8"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$5"\t"$6"\t"$7}' ../data/red.${window}.INTERGENIC.rec.resampling.Z.sex2.CMPAR.txt | paste - ../data/red.${window}.succ.fail > ../data/red.${window}.INTERGENIC.rec.resampling.Z.forR.txt
done

Rscript Ostrich_Z_polymorphism.R 
#________________________________________END CODE BLOCK 10____________________________________#
#___________________________________________SCRIPT NAME_______________________________________#
#___________________________________________R_scripts.sh______________________________________#
#_____________________________________________RUN_____________________________________________#
#_______________________________________bash R_scripts.sh_____________________________________#
#_____________________________________________________________________________________________#
#_____________________________________________________________________________________________#


#!/usr/bin/bash
#______________________________________START CODE BLOCK 11____________________________________#
mkdir -p ../data/female_data
module load bioinfo-tools vcftools BEDTools tabix
echo "Split for females for each subspecies and filter"
zvcf='all.bqsr.gatk.annotated.snp.recal.Z.vcf.gz'    
for subspecies in black_female.txt blue_female.txt red_female.txt # black_female.txt blue_female.txt red_female.txt are located in ostrich_Z_diversity/data/
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
vcftools --gzvcf ../data/${zvcf} --keep ../data/${subspecies} --remove-filtered-all --remove-filtered-geno-all --max-alleles 2 --max-missing 1.0 --maf 0.05 --hwe 0.05 --recode --recode-INFO-all --stdout > ../data/female_data/${vcfname}
bgzip ../data/female_data/${vcfname}
tabix -p vcf ../data/female_data/${vcfname}.gz
doneq

echo "Take allele counts of each subspecies"
for species in black blue red
do
vcftools --gzvcf ../data/female_data/${species}_female.bqsr.gatk.annotated.snp.recal.Z.vcf.gz --counts --out ../data/female_data/${species}.female.counts
done
#*********************************************************************************************#
for subspecies in black blue red
do
zgrep -v "##" ../data/female_data/${subspecies}_female.bqsr.gatk.annotated.snp.recal.Z.vcf.gz | awk '{print $1"\t"$2"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' | awk '{split($3, a, ":") split ($4, b, ":") split ($5, c, ":") split ($6, d, ":") split ($7, e, ":"); print $1"\t"$2"\t"a[1]"\t"b[1]"\t"c[1]"\t"d[1]"\t"e[1]}' | sed 's/\///g' > ../data/female_data/${subspecies}.female.genotypes
paste ../data/female_data/${subspecies}.female.counts.frq.count ../data/female_data/${subspecies}.female.genotypes > ../data/female_data/${subspecies}.female.counts.genotypes 
done

for subspecies in black blue red
do
echo ${subspecies}
python female_nonPAR.py ../data/female_data/${subspecies}.female.counts.genotypes > ../data/female_data/${subspecies}.female.counts.genotypes.nonPAR.recoded
python pi_dxy_female.py ../data/female_data/${subspecies}.female.counts.genotypes.nonPAR.recoded ../data/female_data/${subspecies}.female.counts.success.failure.frq.count > ../data/female_data/${subspecies}.female.pi.txt
done
#*************************************************************************************************************************************************************************************************************************
# Get sum of pi per window for females using male data
for window in 100Kb 200Kb 500Kb 1000Kb
do
echo ${window}
for subspecies in black blue red
do
echo ${subspecies}
python window_female.py ../data/${subspecies}.${window}.repeatCDS.overlap.density.sorted.pi.INTERGENIC.txt ../data/female_data/${subspecies}.female.pi.txt | sort -k2n,2 > ../data/female_data/${subspecies}.${window}.FEMALE.repeatCDS.overlap.density.sorted.pi.INTERGENIC.txt
done
done

# Add confidence interval for intergenic regions
for window in 500Kb 1000Kb
do
echo ${window}
for subspecies in black blue red
do
echo ${subspecies}
python resampling_scaffold.py ../data/female_data/${subspecies}.${window}.FEMALE.repeatCDS.overlap.density.sorted.pi.INTERGENIC.txt  > ../data/female_data/${subspecies}.${window}.female.INTERGENIC.mean.CI.txt
python scaffold_to_chr_vcf.py ../data/female_data/${subspecies}.${window}.female.INTERGENIC.mean.CI.txt > ../data/female_data/${subspecies}.${window}.female.INTERGENIC.mean.CI.ChrZ.txt
done
done

# Fst between males and females
# Create a file with males and females
for subspecies in black blue red
do 
python male_fem_allele_count.py ../data/${subspecies}.male.counts.frq.count ../data/female_data/${subspecies}.${window}.FEMALE.repeatCDS.overlap.density.sorted.pi.INTERGENIC.txt | sort -k1,1 -k2n,2 > ../data/${subspecies}.male.female.counts.Forfst.txt 
done

module load bioinfo-tools vcftools

vcftools --gzvcf black.bqsr.gatk.annotated.snp.recal.Z.vcf.gz --weir-fst-pop black_male.txt --weir-fst-pop black_female.txt --fst-window-size 100000 --out black_male_female_100Kb
vcftools --gzvcf black.bqsr.gatk.annotated.snp.recal.Z.vcf.gz --weir-fst-pop black_male.txt --weir-fst-pop black_female.txt --fst-window-size 200000 --out black_male_female_200Kb

vcftools --gzvcf blue.bqsr.gatk.annotated.snp.recal.Z.vcf.gz --weir-fst-pop blue_male.txt --weir-fst-pop blue_female.txt --fst-window-size 100000 --out blue_male_female_100Kb
vcftools --gzvcf blue.bqsr.gatk.annotated.snp.recal.Z.vcf.gz --weir-fst-pop blue_male.txt --weir-fst-pop blue_female.txt --fst-window-size 200000 --out blue_male_female_200Kb

vcftools --gzvcf red.bqsr.gatk.annotated.snp.recal.Z.vcf.gz --weir-fst-pop red_male.txt --weir-fst-pop red_female.txt --fst-window-size 100000 --out red_male_female_100Kb
vcftools --gzvcf red.bqsr.gatk.annotated.snp.recal.Z.vcf.gz --weir-fst-pop red_male.txt --weir-fst-pop red_female.txt --fst-window-size 200000 --out red_male_female_200Kb

# Convert scaffold to chromosome coordinates
python scaffold_to_chr_vcf.py ../data/black_male_female_100Kb.windowed.weir.fst > ../data/black_male_female_100Kb.windowed.weir.Z.fst
python scaffold_to_chr_vcf.py ../data/black_male_female_200Kb.windowed.weir.fst > ../data/black_male_female_200Kb.windowed.weir.Z.fst

python scaffold_to_chr_vcf.py ../data/blue_male_female_100Kb.windowed.weir.fst > ../data/blue_male_female_100Kb.windowed.weir.Z.fst
python scaffold_to_chr_vcf.py ../data/blue_male_female_200Kb.windowed.weir.fst > ../data/blue_male_female_200Kb.windowed.weir.Z.fst

python scaffold_to_chr_vcf.py ../data/red_male_female_100Kb.windowed.weir.fst > ../data/red_male_female_100Kb.windowed.weir.Z.fst
python scaffold_to_chr_vcf.py ../data/red_male_female_200Kb.windowed.weir.fst > ../data/red_male_female_200Kb.windowed.weir.Z.fst

vcftools --gzvcf all.bqsr.gatk.annotated.snp.recal.Z.vcf.gz --weir-fst-pop black_male.txt --weir-fst-pop blue_male.txt --fst-window-size 200000 --out black_blue_male_200Kb
vcftools --gzvcf all.bqsr.gatk.annotated.snp.recal.Z.vcf.gz --weir-fst-pop black_male.txt --weir-fst-pop red_male.txt --fst-window-size 200000 --out black_red_male_200Kb
vcftools --gzvcf all.bqsr.gatk.annotated.snp.recal.Z.vcf.gz --weir-fst-pop blue_male.txt --weir-fst-pop red_male.txt --fst-window-size 200000 --out blue_red_male_200Kb

python scaffold_to_chr_vcf.py ../data/black_blue_male_200Kb.windowed.weir.fst > ../data/black_blue_male_200Kb.Z
python scaffold_to_chr_vcf.py ../data/black_red_male_200Kb.windowed.weir.fst > ../data/black_red_male_200Kb.Z
python scaffold_to_chr_vcf.py ../data/blue_red_male_200Kb.windowed.weir.fst > ../data/blue_red_male_200Kb.Z
#________________________________________END CODE BLOCK 11____________________________________#
#___________________________________________SCRIPT NAME_______________________________________#
#_______________________________________male_female_fst.sh____________________________________#
#_____________________________________________RUN_____________________________________________#
#_____________________________________bash male_female_fst.sh_________________________________#
#_____________________________________________________________________________________________#
#_____________________________________________________________________________________________#

#!/usr/bin/bash
#______________________________________START CODE BLOCK 12____________________________________#
# Calculate Tajima's D for SNPs that have passed the filtering criteria and are in intergenic
# region using the following steps:

# 1. Use BEDTools to Overlap the vcf with the following file: 
# ${species}.${window}.repeatCDS.overlap.density.sorted.pi.INTERGENIC.txt
# 2. Use vcftools to calculate Tajima's D for 100 Kb and 200 Kb windows
# 3. Assign scaffolds to chromosomes

# Load the necessary modules
module load bioinfo-tools vcftools BEDTools tabix

for species in black blue red
do
echo ${species}
#echo "intersect vcf and intergenic"
#awk '{print $1"\t"$2"\t"$3}' ../data/${species}.100Kb.repeatCDS.overlap.density.sorted.pi.INTERGENIC.txt > ../data/${species}.intergenic.forvcf.bed
#cat header ../data/${species}.intergenic.forvcf.bed > ../data/${species}.intergenic.forvcf.header.bed
#vcftools --gzvcf ../data/${species}_male.bqsr.gatk.annotated.snp.recal.Z.vcf.gz --bed ../data/${species}.intergenic.forvcf.header.bed --recode --recode-INFO-all --stdout > ../data/${species}_male.INTERGENIC.recal.Z.vcf

#echo "bgzip and tabix"
#bgzip ../data/${species}_male.INTERGENIC.recal.Z.vcf
#tabix -p vcf ../data/${species}_male.INTERGENIC.recal.Z.vcf.gz

echo "get Tajima's D"
vcftools --gzvcf ../data/${species}_male.bqsr.gatk.annotated.snp.recal.Z.vcf.gz --TajimaD 100000 --out ../data/${species}_male.100Kb.allsites
vcftools --gzvcf ../data/${species}_male.bqsr.gatk.annotated.snp.recal.Z.vcf.gz --TajimaD 200000 --out ../data/${species}_male.200Kb.allsites

python scaffold_to_chr_vcf.py ../data/${species}_male.100Kb.allsites.Tajima.D > ../data/${species}_male_100Kb.allsites.Tajima.D.chrZ.txt
python scaffold_to_chr_vcf.py ../data/${species}_male.200Kb.allsites.Tajima.D > ../data/${species}_male_200Kb.allsites.Tajima.D.chrZ.txt
done
#________________________________________END CODE BLOCK 12____________________________________#
#___________________________________________SCRIPT NAME_______________________________________#
#___________________________________________tajima_D.sh_______________________________________#
#______________________________________________RUN____________________________________________#
#________________________________________bash tajima_D.sh_____________________________________#
#_____________________________________________________________________________________________#
#_____________________________________________________________________________________________#