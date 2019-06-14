#!/usr/bin/python
from __future__ import division
import sys

f1 = open(sys.argv[1], "r")

scaffold_order = ["superscaffold26", "superscaffold54", "superscaffold35", "superscaffold36", "superscaffold62", "superscaffold67", "superscaffold69-1", "superscaffold93", "superscaffold63", "superscaffold88", "superscaffold83", "superscaffold92"]

f1_dict = {}
f1_list = []
for line in f1:
	line = line.strip("\n").split("\t")
	key = line[0]+"_"+line[12]+"_"+line[13]
	value = line
	if not key in f1_list:
		f1_list.append(key)
	if key in f1_dict.keys():
		f1_dict[key].append(value)
	else:
		f1_dict[key] = [value]

#175224  77196   161059  161087  48809   48868   6147    2657    60070   66414   100986  46979   187155  187157  68091   27560
#print(f1_list)
#print(f1_dict)
header=["Scaffold", "Window_start", "Window_end", "pi_per_window", "GC_density", "repeat_density", "CDS_density"]
print("\t".join(header))
for scaffold in scaffold_order:
	for key in f1_list:
		if key.split("_")[0] == scaffold:
			#print(key)
		#for key in f1_dict.keys():
			# Number of segregating sites per window
			#snp_count = len(f1_dict[key])
			pi_count = []
			for element in f1_dict[key]:
				pi_count.append(float(element[6]))
			pi_count_key = sum(pi_count)
			if key.split("_")[0] == "superscaffold54":
				if int(key.split("_")[1]) > 16379243: # Removing the non-Z linked part of superscaffold54
					pass
				else:
					print(key.split("_")[0]+"\t"+key.split("_")[1]+"\t"+key.split("_")[2]+"\t"+str(pi_count_key/int(element[14]))+"\t"+str(int(element[15])/int(element[14]))+"\t"+str((int(element[17])-int(element[16]))/int(element[14]))+"\t"+str(int(element[22])/int(element[14])))  # repeat density# CDS density
			else:
			#print(element[11])
				print(key.split("_")[0]+"\t"+key.split("_")[1]+"\t"+key.split("_")[2]+"\t"+str(pi_count_key/int(element[14]))+"\t"+str(int(element[15])/int(element[14]))+"\t"+str((int(element[17])-int(element[16]))/int(element[14]))+"\t"+str(int(element[22])/int(element[14])))  # repeat density# CDS density


