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