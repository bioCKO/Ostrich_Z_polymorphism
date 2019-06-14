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
cut -f4 ../data/${subspecies}_rec_female_per_${window}.txt | awk 'BEGIN{print "female_CM_per_bp"}{if(NR>1) print $1}' | paste ../data/${subspecies}.${window}.INTERGENIC.rec.resampling.Z.sex.txt - > ../data/${subspecies}.${window}.INTERGENIC.rec.resampling.Z.sex2.txt   
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