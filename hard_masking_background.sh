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