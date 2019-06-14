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