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