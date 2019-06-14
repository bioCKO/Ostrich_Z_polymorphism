#!/usr/bin/python
import sys


# Written by Homa Papoli - 22 Nov 2018
# There is a problem in formatting of optical mapping gff of
# ostrich. Genes in the reverse strand should be sorted from 
# the largest coordinate to the smallest as of the format of 
# other gffs. There are some inconsistency in case of ostrich which
# can get corrected with the help of the following script. 

gff = open("../data/Struthio_camelus.OM.gene.20130116.gff", "r")

mRNA_dict = {}
gff_dict = {}
for line in gff:
	line = line.strip("\n").split("\t")
	#print(line[0])
	if line[2] == "CDS":
		key = line[0]+":"+line[8].split("=")[1].replace(";", "")
		value = [line[3], line[4], line[6], line[7], line[8], line[0], line[1], line[2], line[5]]
		value = [line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8]]
		if key in gff_dict.keys():
			gff_dict[key].append(value)
		else:
			gff_dict[key] = [value]
	elif line[2] == "mRNA":
		key = line[0]+":"+line[8].split(";")[0].split("=")[1]
		value = [line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8]]
		if key in mRNA_dict.keys():
			mRNA_dict[key].append(value)
		else:
			mRNA_dict[key] = [value]		

#print(gff_dict)
#print(mRNA_dict)

for key in mRNA_dict.keys():
	print("\t".join(mRNA_dict[key][0]))
#for key in gff_dict.keys():
	start_pos = [int(i[3]) for i in gff_dict[key]]
	end_pos = [int(i[4]) for i in gff_dict[key]]
	if gff_dict[key][0][6] == "+":
		if not start_pos == sorted(start_pos):
			#print(key, gff_dict[key][0][2], start_pos, end_pos)
			#print(sorted(start_pos))
			for index, item in enumerate(sorted(start_pos)):
				l = [gff_dict[key][index][0], gff_dict[key][index][1], gff_dict[key][index][2], item, sorted(end_pos)[index], gff_dict[key][index][5], gff_dict[key][index][6], gff_dict[key][index][7], gff_dict[key][index][8]]
				print("\t".join(str(i) for i in l))
		else:
			for index, item in enumerate(start_pos):
				l = [gff_dict[key][index][0], gff_dict[key][index][1], gff_dict[key][index][2], item, end_pos[index], gff_dict[key][index][5], gff_dict[key][index][6], gff_dict[key][index][7], gff_dict[key][index][8]]
				print("\t".join(str(i) for i in l))
	elif gff_dict[key][0][6] == "-":
		if not start_pos == sorted(start_pos , reverse=True):
			#print(key, gff_dict[key][0][2], start_pos, end_pos)
			for index, item in enumerate(sorted(start_pos , reverse=True)):
				l = [gff_dict[key][index][0], gff_dict[key][index][1], gff_dict[key][index][2], item, sorted(end_pos, reverse=True)[index], gff_dict[key][index][5], gff_dict[key][index][6], gff_dict[key][index][7], gff_dict[key][index][8]]
				print("\t".join(str(i) for i in l))
		else:
			for index, item in enumerate(start_pos):
				l = [gff_dict[key][index][0], gff_dict[key][index][1], gff_dict[key][index][2], item, end_pos[index], gff_dict[key][index][5], gff_dict[key][index][6], gff_dict[key][index][7], gff_dict[key][index][8]]
				print("\t".join(str(i) for i in l))

# Print out the reformated gff
