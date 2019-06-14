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