#!/usr/bin/python
import sys
from collections import defaultdict

#************************************************************************************
# Written by Homa Papoli 
# The script converts the scaffold based coordinates of postion on ostrich Z to 
# Z chromosome coordinates based on linkage map obtained in Yazdi and Ellegren 2018. 
#************************************************************************************

# run ./scaffold_to_chr_vcf.py ${dirs}/black.bqsr.gatk.annotated.snp.recal.vcf.Z.thinned1K > toto2

infile = open(sys.argv[1], "r")

scaffold_order = ["superscaffold26", "superscaffold54", "superscaffold35", "superscaffold36", "superscaffold62", "superscaffold67", "superscaffold69-1", "superscaffold93", "superscaffold63", "superscaffold88", "superscaffold83", "superscaffold92"]

# superscaffold54 length pre LM: 29256470
scaffold_length = {'superscaffold26': 25310599, 'superscaffold54':16379243 , 'superscaffold35': 4625539, 'superscaffold36': 9394175, 'superscaffold62': 2917291, 'superscaffold67': 5300260, 'superscaffold69-1': 5978518, 'superscaffold93': 4983591, 'superscaffold63': 1692925, 'superscaffold88': 624114, 'superscaffold83': 782506, 'superscaffold92': 2882843}

# Read vcf file into a dictionary
header = []
vcf_dict = {}
for line in infile:
	if line.startswith("Scaffold") or line.startswith("CHROM"):
		header.append(line.strip("\n"))
	elif not line.startswith("#"):
		line = line.strip("\n").split()
		key, value = line[0], line[1:]
		if key in vcf_dict.keys():
			vcf_dict[key].append(value)
		else:
			vcf_dict[key] = [value]
#print(vcf_dict)			
# Print header, adding ChrZ
#for l in header:
	#if not l.startswith("##reference"):
	#	if not l.startswith("##source"):
	#		if not l.startswith("#CHROM"):
	#print l
#print "##contig=<ID=ChrZ,length=80871604>"		
#for l in header:
#	if l.startswith("##reference") or l.startswith("##source") or l.startswith("#CHROM"):
#		print l
		
#print("\t".join(header))
# Change the coordinates
for scaffold in scaffold_order:
	if scaffold == "superscaffold26":
		for l in vcf_dict[scaffold]:
			print ("ChrZ" + "\t" + scaffold + "\t" + "\t".join(l[0:2]) + "\t" + "\t".join(l))
# ******************************************************************************			
	elif scaffold == "superscaffold54":
		new_list = []
		scaf_list = []
		for l in vcf_dict[scaffold]:
			if int(l[0]) <= scaffold_length[scaffold]:
				scaf_list.append(l[0:2])
				l[0] = str(scaffold_length["superscaffold26"] + 5000 + (scaffold_length[scaffold] - int(l[0]) + 1))
				l[1] = str(scaffold_length["superscaffold26"] + 5000 + (scaffold_length[scaffold] - int(l[1]) + 1))
				new_list.append(l)
		for i in reversed(new_list):
			ind = new_list.index(i)
			print ("ChrZ" + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(i))
# ******************************************************************************			
	elif scaffold == "superscaffold35":
		scaf_list = []
		for l in vcf_dict[scaffold]:
			ind = vcf_dict[scaffold].index(l)
			scaf_list.append(l[0:2])
			l[0] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + int(l[0]))
			l[1] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + int(l[1]))
			print ("ChrZ" + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l))
# ******************************************************************************			
	elif scaffold == "superscaffold36":
		new_list = []
		scaf_list = []
		for l in vcf_dict[scaffold]:
			scaf_list.append(l[0:2])
			l[0] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + (scaffold_length[scaffold] - int(l[0]) + 1))
			l[1] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + (scaffold_length[scaffold] - int(l[1]) + 1))
			new_list.append(l)
		for i in reversed(new_list):
			ind = new_list.index(i)
			print ("ChrZ" + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(i))
# ******************************************************************************			
	elif scaffold == "superscaffold62":
		scaf_list = []
		for l in vcf_dict[scaffold]:
			ind = vcf_dict[scaffold].index(l)
			scaf_list.append(l[0:2])
			l[0] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + int(l[0]))
			l[1] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + int(l[1]))
			print ("ChrZ" + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t"+ "\t".join(l))
# ******************************************************************************			
	elif scaffold == "superscaffold67":
		scaf_list = []
		for l in vcf_dict[scaffold]:
			ind = vcf_dict[scaffold].index(l)
			scaf_list.append(l[0:2])
			l[0] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + scaffold_length["superscaffold62"] + 5000 + int(l[0]))
			l[1] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + scaffold_length["superscaffold62"] + 5000 + int(l[1]))			
			print ("ChrZ" + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l))
# ******************************************************************************			
	elif scaffold == "superscaffold69-1":
		scaf_list = []
		for l in vcf_dict[scaffold]:
			ind = vcf_dict[scaffold].index(l)
			scaf_list.append(l[0:2])
			l[0] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + scaffold_length["superscaffold62"] + 5000 + scaffold_length["superscaffold67"] + 5000 + int(l[0]))
			l[1] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + scaffold_length["superscaffold62"] + 5000 + scaffold_length["superscaffold67"] + 5000 + int(l[1]))			
			print ("ChrZ" + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l))
# ******************************************************************************			
	elif scaffold == "superscaffold93":
		scaf_list = []
		for l in vcf_dict[scaffold]:
			ind = vcf_dict[scaffold].index(l)
			scaf_list.append(l[0:2])
			l[0] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + scaffold_length["superscaffold62"] + 5000 + scaffold_length["superscaffold67"] + 5000 + scaffold_length["superscaffold69-1"] + 5000 + int(l[0]))
			l[1] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + scaffold_length["superscaffold62"] + 5000 + scaffold_length["superscaffold67"] + 5000 + scaffold_length["superscaffold69-1"] + 5000 + int(l[1]))			
			print ("ChrZ" + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l))
# ******************************************************************************			
	elif scaffold == "superscaffold63":
		scaf_list = []
		for l in vcf_dict[scaffold]:
			ind = vcf_dict[scaffold].index(l)
			scaf_list.append(l[0:2])
			l[0] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + scaffold_length["superscaffold62"] + 5000 + scaffold_length["superscaffold67"] + 5000 + scaffold_length["superscaffold69-1"] + 5000 + scaffold_length["superscaffold93"] + 5000 + int(l[0]))
			l[1] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + scaffold_length["superscaffold62"] + 5000 + scaffold_length["superscaffold67"] + 5000 + scaffold_length["superscaffold69-1"] + 5000 + scaffold_length["superscaffold93"] + 5000 + int(l[1]))			
			print ("ChrZ" + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l))
# ******************************************************************************			
	elif scaffold == "superscaffold88":
		scaf_list = []
		for l in vcf_dict[scaffold]:
			ind = vcf_dict[scaffold].index(l)
			scaf_list.append(l[0:2])
			l[0] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + scaffold_length["superscaffold62"] + 5000 + scaffold_length["superscaffold67"] + 5000 + scaffold_length["superscaffold69-1"] + 5000 + scaffold_length["superscaffold93"] + 5000 + scaffold_length["superscaffold63"] + 5000 + int(l[0]))
			l[1] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + scaffold_length["superscaffold62"] + 5000 + scaffold_length["superscaffold67"] + 5000 + scaffold_length["superscaffold69-1"] + 5000 + scaffold_length["superscaffold93"] + 5000 + scaffold_length["superscaffold63"] + 5000 + int(l[1]))
			print ("ChrZ" + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l))
# ******************************************************************************			
	elif scaffold == "superscaffold83":
		scaf_list = []
		for l in vcf_dict[scaffold]:
			ind = vcf_dict[scaffold].index(l)
			scaf_list.append(l[0:2])
			l[0] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + scaffold_length["superscaffold62"] + 5000 + scaffold_length["superscaffold67"] + 5000 + scaffold_length["superscaffold69-1"] + 5000 + scaffold_length["superscaffold93"] + 5000 + scaffold_length["superscaffold63"] + 5000 + scaffold_length["superscaffold88"] + 5000 + int(l[0]))
			l[1] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + scaffold_length["superscaffold62"] + 5000 + scaffold_length["superscaffold67"] + 5000 + scaffold_length["superscaffold69-1"] + 5000 + scaffold_length["superscaffold93"] + 5000 + scaffold_length["superscaffold63"] + 5000 + scaffold_length["superscaffold88"] + 5000 + int(l[1]))			
			print ("ChrZ" + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l))
# ******************************************************************************			
	elif scaffold == "superscaffold92":
		scaf_list = []
		for l in vcf_dict[scaffold]:
			ind = vcf_dict[scaffold].index(l)
			scaf_list.append(l[0:2])
			l[0] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + scaffold_length["superscaffold62"] + 5000 + scaffold_length["superscaffold67"] + 5000 + scaffold_length["superscaffold69-1"] + 5000 + scaffold_length["superscaffold93"] + 5000 + scaffold_length["superscaffold63"] + 5000 + scaffold_length["superscaffold88"] + 5000 + scaffold_length["superscaffold83"] + 5000 + int(l[0]))
			l[1] = str(scaffold_length["superscaffold26"] + 5000 + scaffold_length["superscaffold54"] + 5000 + scaffold_length["superscaffold35"] + 5000 + scaffold_length["superscaffold36"] + 5000 + scaffold_length["superscaffold62"] + 5000 + scaffold_length["superscaffold67"] + 5000 + scaffold_length["superscaffold69-1"] + 5000 + scaffold_length["superscaffold93"] + 5000 + scaffold_length["superscaffold63"] + 5000 + scaffold_length["superscaffold88"] + 5000 + scaffold_length["superscaffold83"] + 5000 + int(l[1]))			
			print ("ChrZ" + "\t" + scaffold + "\t" + ["\t".join(map(str,x)) for x in scaf_list][ind] + "\t" + "\t".join(l))
# ******************************************************************************				
	
				        
