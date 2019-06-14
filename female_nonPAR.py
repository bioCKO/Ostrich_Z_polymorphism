#!/usr/bin/python
import sys

# To run:
# python female_Z_nonPAR_het.py ../data/female_data/black.female.counts.genotypes ../data/female_data/nonPAR_outfemale ../data/female_data/PAR_outfemale ../data/female_data/scaf36_outfemale ../data/female_data/scaf62_outfemale


f1 = open(sys.argv[1], "r")

scaffold_order = ["superscaffold26", "superscaffold54", "superscaffold35", "superscaffold36", "superscaffold62", "superscaffold67", "superscaffold69-1", "superscaffold93", "superscaffold63", "superscaffold88", "superscaffold83", "superscaffold92"]
nonPAR = ["superscaffold36", "superscaffold62", "superscaffold67", "superscaffold69-1", "superscaffold93", "superscaffold63", "superscaffold88", "superscaffold83", "superscaffold92"]

for line in f1:
	line = line.strip("\n").split("\t")
	if not line[0] in nonPAR: # If the scaffold is in the recombining section of Z in females:
		print("\t".join(line[0:6])) # Print as it is
	else: # If the scaffold is in the non-recombining section of Z in females:
		if not line[0] == "superscaffold36":
			allele1 = line[4].split(":")[1]
			if int(allele1)%2 == 0:
				if len(line) == 12:
					if line[8]!="01" and line[9]!="01" and line[10]!="01" and line[11]!="01":
						print(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+str(int(line[3])//2)+"\t"+line[4].split(":")[0]+":"+str(int(line[4].split(":")[1])//2)+"\t"+line[5].split(":")[0]+":"+str(int(line[5].split(":")[1])//2))
				elif len(line) == 13:
					if line[8]!="01" and line[9]!="01" and line[10]!="01" and line[11]!="01" and line[12]!="01":
						print(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+str(int(line[3])//2)+"\t"+line[4].split(":")[0]+":"+str(int(line[4].split(":")[1])//2)+"\t"+line[5].split(":")[0]+":"+str(int(line[5].split(":")[1])//2))
				else:
					raiseException("Line columns don't match expectation")
			else:
				pass
		if line[0] == "superscaffold36":
			if int(line[1]) > 5408561:
				print("\t".join(line[0:6]))
			else:
				allele1 = line[4].split(":")[1]
				if int(allele1)%2 == 0:
					if len(line) == 12:
						if line[8]!="01" and line[9]!="01" and line[10]!="01" and line[11]!="01":
							print(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+str(int(line[3])//2)+"\t"+line[4].split(":")[0]+":"+str(int(line[4].split(":")[1])//2)+"\t"+line[5].split(":")[0]+":"+str(int(line[5].split(":")[1])//2))
					elif len(line) == 13:
						if line[8]!="01" and line[9]!="01" and line[10]!="01" and line[11]!="01" and line[12]!="01":
							print(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+str(int(line[3])//2)+"\t"+line[4].split(":")[0]+":"+str(int(line[4].split(":")[1])//2)+"\t"+line[5].split(":")[0]+":"+str(int(line[5].split(":")[1])//2))
					else:
						raiseException("Line columns don't match expectation")					
