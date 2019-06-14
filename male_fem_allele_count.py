#!/usr/bin/python
import sys

male = open(sys.argv[1], "r")
female = open(sys.argv[2], "r")

scaffold_order = ["superscaffold26", "superscaffold54", "superscaffold35", "superscaffold36", "superscaffold62", "superscaffold67", "superscaffold69-1", "superscaffold93", "superscaffold63", "superscaffold88", "superscaffold83", "superscaffold92"]

male_dict = {}
for line in male:
	line = line.strip("\n").split("\t")
	if not line[0] == "CHROM":
		key = line[0]+"_"+line[1]
		value = line[4:]
		male_dict[key] = value
#print(male_dict)
female_dict = {}
for line in female:
	line = line.strip("\n").split("\t")
	if not line[0] == "CHROM":
		key = line[0]+"_"+line[2]
		value = line[4:6]
		female_dict[key] = value
#print(female_dict)
for scaffold in scaffold_order:
	for key in female_dict.keys():
		if key.split("_")[0] == scaffold:
			if key in male_dict.keys():
				print(key.split("_")[0]+"\t"+key.split("_")[1]+"\t"+male_dict[key][0].split(":")[1]+"\t"+male_dict[key][1].split(":")[1]+"\t"+female_dict[key][0].split(":")[1]+"\t"+female_dict[key][1].split(":")[1])	
