#!/usr/bin/python
from __future__ import division
import sys

#***********************************************************************
# Written by Homa Papoli - Nov 2018
#***********************************************************************
# Calculate the rate of recombination by obtaining the slope of the 
# graph, genetic position as a function of physical position. This
# considers a linear relationship between genetic position and physical
# position. While this is violated at large physical distances, at a 
# short distance, the slope of the curve can be approximated by a line.
#***********************************************************************
# Usage: ./recombination.py ../data/LGZ3.sex_averaged.lifted.bed
#***********************************************************************

pos_l = []
rec_f = open(sys.argv[1], "r")
for line in rec_f:
	line = line.strip("\n").split("\t")
#	print(line)
	element = line[4]+":"+line[3].split(":")[1]+":"+line[1]+":"+line[2]
	#print(element)
	pos_l.append(element)

for index, item in enumerate(pos_l):
	if index == 0:
		print(item.split(":")[0]+"\t"+"0"+"\t"+item.split(":")[1]+"\t"+item.split(":")[2]+"\t"+"0"+"\t"+pos_l[index].split(":")[4])
	else:
		if pos_l[index].split(":")[0] == pos_l[index-1].split(":")[0]: # Calculate the recombinate rate only within scaffolds.
			rec_rate = abs(float(pos_l[index-1].split(":")[2])-float(pos_l[index].split(":")[2]))/abs(int(pos_l[index-1].split(":")[1])-int(pos_l[index].split(":")[1]))
			print(item.split(":")[0]+"\t"+pos_l[index-1].split(":")[1]+"\t"+pos_l[index].split(":")[1]+"\t"+str(rec_rate)+"\t"+pos_l[index-1].split(":")[4]+"\t"+pos_l[index].split(":")[4])
		else:
			pass
