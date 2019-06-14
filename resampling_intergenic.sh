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