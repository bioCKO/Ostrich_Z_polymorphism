#!/usr/bin/python
from __future__ import division
import sys

#****************************************************************************
# Written by Homa Papoli - Jan 2019
#****************************************************************************
# This script prints window start and end for chromosomal coordinates and the
# midpoint. 
# Run:
# python CM_to_PAR_boundary.py  ../data/black.${window}.INTERGENIC.rec.resampling.Z.sex2.txt > ../data/black.${window}.INTERGENIC.rec.resampling.Z.sex4.txt
#****************************************************************************

f1 = open(sys.argv[1], "r")

PAR_boundary = 53065700
f1_dict = {}
for line in f1:
	line = line.strip("\n").split("\t")
	if not line[0] == "Scaffold":
		chrZ_start = int(line[8])
		chrZ_end = int(line[9])
		midpoint = (chrZ_start + chrZ_end)/2
		if midpoint < PAR_boundary:
			print("ChrZ"+"\t"+str(int(midpoint))+"\t"+str(PAR_boundary))
		else:
			print("ChrZ"+"\t"+str(PAR_boundary)+"\t"+str(int(midpoint)))