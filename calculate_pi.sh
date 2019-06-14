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