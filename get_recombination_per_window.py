#!/usr/bin/python
from __future__ import division
import sys

#****************************************************************************
# Written by Homa Papoli - Jan 2019
#****************************************************************************
# Script calculates recombination rate per window in cM per base pair (cM/bp)
# Logic:
# Given recombination rate per site between Marker 1 and Marker 2:
# _________|_______________|__________
#       Marker 1        Marker 2
# Case 1:
# _________|_______________|__________
#       Marker 1        Marker 2
#            |__________|
#	 Window_start    Window_end
# Real example:
# 'superscaffold62:1000000:2000000': [['superscaffold62', '207526', '2434793', '0.0', '55917482', '58144749', '1000000']]
# Calculation: 1000000*0.0/1000000
# Case 2:
#				r1		  r2
# _________|_________|_________|_
#       Marker 1   Marker 2  Marker 3
#            |__________|
#	 Window_start    Window_end
# Calculation: r1*(Overlap of window with the distance between Marker 1 and Marker 2) + r2*(Overlap of window with the distance between Marker 2 and Marker 3)
# Real example:
# 'superscaffold26:3000000:4000000': [['superscaffold26', '1113307', '3461665', '2.19131835947e-06', '1113307', '3461665', '461665'], 
#                                     ['superscaffold26', '3461665', '4790876', '4.25891750821e-06', '3461665', '4790876', '538335']]
# Calculation: (461665*2.19131835947e-06) + (538335*4.25891750821e-06) = 3.304379347206948
# Recombination rate: 3.304379347206948/(461665+538335) = 3.304379347206948e-06
#*****************************************************************************

f1 = open(sys.argv[1], "r")
f2 = open(sys.argv[2], "r")

f2_l = []
for line in f2:
	if not line.startswith("Scaffold"):
		line = line.strip("\n").split("\t")
		f2_l.append(line[0]+":"+line[1]+":"+line[2])

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
print("Scaffold"+"\t"+"Window_start"+"\t"+"Window_end"+"\t"+"CM_per_bp")
for key in f2_l:
	if len(f1_dict[key]) == 1:
		if f1_dict[key][0][0] == ".":
			print(key.split(":")[0]+"\t"+key.split(":")[1]+"\t"+key.split(":")[2]+"\t"+"NA")
		else:
			print(key.split(":")[0]+"\t"+key.split(":")[1]+"\t"+key.split(":")[2]+"\t"+f1_dict[key][0][3])
	else:
		addup = 0
		total = 0 #int(key.split(":")[2])-int(key.split(":")[1])
		for element in f1_dict[key]:
			addup = addup + int(element[6])*float(element[3])
			total = total + int(element[6])
		new_rho = addup/total
		print(key.split(":")[0]+"\t"+key.split(":")[1]+"\t"+key.split(":")[2]+"\t"+str(new_rho))






