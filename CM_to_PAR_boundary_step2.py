#!/usr/bin/python
from __future__ import division
import sys

#****************************************************************************
# Written by Homa Papoli - Jan 2019
#****************************************************************************
# This scripts calculates the genetic distance of the focal point, i.e. the
# mid section of a window from the PAR boundary by summing of the genetic 
# distances between the markers on the chromosome. 
# Run: python CM_to_PAR_boundary_step2.py ../data/black.${window}.INTERGENIC.rec.overlap.txt test4 > ../data/distance_from_PAR_${window}.txt
#****************************************************************************

f1 = open(sys.argv[1], "r")
f2 = open(sys.argv[2], "r")

f2_l = []
for line in f2:
	if not line.startswith("Scaffold"):
		line = line.strip("\n").split("\t")
		f2_l.append(line[0]+":"+line[1]+":"+line[2])
#print(f2_l)

f1_dict = {}
for line in f1:
	line = line.strip("\n").split("\t")
	key = line[0]+":"+line[1]+":"+line[2]
	value = line[3:]
	if key in f1_dict.keys():
		f1_dict[key].append(value)
	else:
		f1_dict[key] = [value]

#print(f1_dict)

print("Scaffold"+"\t"+"Focal_site"+"\t"+"PAR_boundary"+"\t"+"CM_from_PAR_boundary")
for key in f1_dict.keys():
	if len(f1_dict[key]) == 1:
		if f1_dict[key][0][0] == ".":
			print(key.split(":")[0]+"\t"+key.split(":")[1]+"\t"+key.split(":")[2]+"\t"+"NA")
		else:
			print(key.split(":")[0]+"\t"+key.split(":")[1]+"\t"+key.split(":")[2]+"\t"+f1_dict[key][0][3])
	else:
		addup = 0
		#total = 0 #int(key.split(":")[2])-int(key.split(":")[1])
		for element in f1_dict[key]:
			if element[3] != "NA":
				addup = addup + int(element[4])*float(element[3])
		#	total = total + int(element[4])
		#new_rho = addup/total
		print(key.split(":")[0]+"\t"+key.split(":")[1]+"\t"+key.split(":")[2]+"\t"+str(addup))






