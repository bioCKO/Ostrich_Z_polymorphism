#!/usr/bin/bash
#______________________________________START CODE BLOCK 2_____________________________________#
# CODE BLOCK 2 deals with:
# Background filtering
# Background DNA is filtered for repeats and coverage.
# |ATG[C/C]GTGTGGNNNNNNNNATGGTGC[A/A]CGTGACGTCTNNNNNNNNNNNNNNNNTGTACACGTG[C/C]GTGACACAACGT|
#                 Repeat                        Coverage filter  
#_____________________________________________________________________________________________#
# Load the necessary modules
module load bioinfo-tools BEDTools
#*********************************************************************************************#
# Get windows per scaffolds
#*********************************************************************************************#
for window in 100Kb 200Kb 500Kb 1000Kb
do
w_size=$(echo $window | sed 's/Kb/000/g')
python sliding_window_ZA.py ../data/Z_scaffolds_length.txt ${w_size} > ../data/ostrich.Z.${window}.bed
done
#*********************************************************************************************#
# Repeat elements
#*********************************************************************************************#
repeat_addr='../data/Ostrich_repeatMask/Struthio_camelus.20130116.OM.fa.out'
awk -v OFS="\t" '$1=$1' ${repeat_addr} | awk '{$6=$6-1; print $5"\t"$6"\t"$7}' | awk 'NR>2' | mergeBed -i - | grep -f ../data/Z_scaffolds.txt > ../data/SC.allRepeats.bed
#*********************************************************************************************#
# Get overlap between windows and element
#*********************************************************************************************#
# Repeats
#*********************************************************************************************#
for window in 100Kb 200Kb 500Kb 1000Kb
do
echo ${window}
bedtools intersect -a ../data/ostrich.Z.${window}.bed -b ../data/SC.allRepeats.bed -wao > ../data/ostrich.${window}.AllRepeat.overlap.txt
done
#*********************************************************************************************#
# Background coverage
#*********************************************************************************************#
coverageAddr='/proj/uppstore2017180/private/cornwallis/data/interim/individual/genomecov'
for window in 100Kb 200Kb 500Kb 1000Kb
do
echo ${window}
for subspecies in black blue red
do
echo ${subspecies}
bedtools intersect -a ../data/ostrich.Z.${window}.bed -b ${coverageAddr}/${subspecies}.coverage.7X.70pct.bed.gz -wao > ../data/${subspecies}.${window}.cov.overlap.txt
done
done
#*********************************************************************************************#
# Sum of overlapping features
#*********************************************************************************************#
for window in 100Kb 200Kb 500Kb 1000Kb
do
echo ${window}
python get_sum_overlapping.py ../data/ostrich.${window}.AllRepeat.overlap.txt > ../data/ostrich.${window}.AllRepeat.overlap.sum.txt
# Use an already repeat masked genome
# Because the repeat masked genome is used, sum of bases and gc for repeat elements is zero in the output
# as due to hard masking, these positions are set to "N" in the repeat masked fasta file. 
python get_density_feature.py ../data/ostrich.${window}.AllRepeat.overlap.sum.txt ../data/Ostrich_repeatMask/Struthio_camelus.20130116.OM.fa.masked ../data/ostrich.${window}.AllRepeat.overlap.density.txt 
done
#*********************************************************************************************#
for window in 100Kb 200Kb 500Kb 1000Kb
do
echo ${window}
for subspecies in black blue red
do
python get_sum_overlapping.py ../data/${subspecies}.${window}.cov.overlap.txt > ../data/${subspecies}.${window}.cov.sum.overlap.txt
# Use an already repeat masked genome
# Because the repeat masked genome is used, sum of bases and gc for repeat elements is zero in the output
# as due to hard masking, these positions are set to "N" in the repeat masked fasta file. 
python get_density_feature.py ../data/${subspecies}.${window}.cov.sum.overlap.txt ../data/Ostrich_repeatMask/Struthio_camelus.20130116.OM.fa.masked ../data/${subspecies}.${window}.cov.overlap.density.txt 
done
done
#________________________________________END CODE BLOCK 2_____________________________________#
#___________________________________________SCRIPT NAME_______________________________________#
#____________________________________background_repeat_cov.sh_________________________________#
#_____________________________________________RUN_____________________________________________#
#________________________________bash background_repeat_cov.sh________________________________#
#_____________________________________________________________________________________________#