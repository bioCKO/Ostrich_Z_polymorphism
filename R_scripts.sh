#!/usr/bin/bash
#______________________________________START CODE BLOCK 10____________________________________#
# Print the columns with an easier-to-read order for statistical analysis
for window in 200Kb 500Kb 1000Kb
do
echo ${window}
awk '{if(NR>1 && $9>$10) print $1"\t"$2"\t"$3"\t"$10"\t"$9"\t"$4"\t"$11"\t"$12"\t"$13"\t"$14"\t"$8"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$5"\t"$6"\t"$7; else print $1"\t"$2"\t"$3"\t"$9"\t"$10"\t"$4"\t"$11"\t"$12"\t"$13"\t"$14"\t"$8"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$5"\t"$6"\t"$7}' ../data/black.${window}.INTERGENIC.rec.resampling.Z.sex2.CMPAR.txt > ../data/black.${window}.INTERGENIC.rec.resampling.Z.forR.txt
awk '{if(NR>1 && $9>$10) print $1"\t"$2"\t"$3"\t"$10"\t"$9"\t"$4"\t"$11"\t"$12"\t"$13"\t"$14"\t"$8"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$5"\t"$6"\t"$7; else print $1"\t"$2"\t"$3"\t"$9"\t"$10"\t"$4"\t"$11"\t"$12"\t"$13"\t"$14"\t"$8"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$5"\t"$6"\t"$7}' ../data/blue.${window}.INTERGENIC.rec.resampling.Z.sex2.CMPAR.txt > ../data/blue.${window}.INTERGENIC.rec.resampling.Z.forR.txt
awk '{if(NR>1 && $9>$10) print $1"\t"$2"\t"$3"\t"$10"\t"$9"\t"$4"\t"$11"\t"$12"\t"$13"\t"$14"\t"$8"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$5"\t"$6"\t"$7; else print $1"\t"$2"\t"$3"\t"$9"\t"$10"\t"$4"\t"$11"\t"$12"\t"$13"\t"$14"\t"$8"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$5"\t"$6"\t"$7}' ../data/red.${window}.INTERGENIC.rec.resampling.Z.sex2.CMPAR.txt > ../data/red.${window}.INTERGENIC.rec.resampling.Z.forR.txt
done

Rscript Ostrich_Z_polymorphism.R 
#________________________________________END CODE BLOCK 8_____________________________________#
#___________________________________________SCRIPT NAME_______________________________________#
#___________________________________________R_scripts.sh______________________________________#
#_____________________________________________RUN_____________________________________________#
#_______________________________________bash R_scripts.sh_____________________________________#
#_____________________________________________________________________________________________#
#_____________________________________________________________________________________________#