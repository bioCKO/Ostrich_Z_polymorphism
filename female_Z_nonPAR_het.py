#!/usr/bin/python
import sys

f1 = open(sys.argv[1], "r")
non_par_out = open(sys.argv[2], "w")
par_out = open(sys.argv[3], "w")
scaf36_out = open(sys.argv[4], "w")
scaf62_out = open(sys.argv[5], "w")

nonPAR = ["superscaffold36", "superscaffold62", "superscaffold67", "superscaffold69-1", "superscaffold93", "superscaffold63", "superscaffold88", "superscaffold83", "superscaffold92"]

# python female_Z_nonPAR_het.py ../data/female_data/black.female.counts.genotypes ../data/female_data/nonPAR_outfemale 
# ../data/female_data/PAR_outfemale ../data/female_data/scaf36_outfemale ../data/female_data/scaf62_outfemale

for line in f1:
	line = line.strip("\n").split("\t")
	if not line[0] == "CHROM":
		if line[0] in nonPAR and line[0] != "superscaffold36" and line[0] != "superscaffold62": # If the scaffold is in the recombining section of Z in females
			genotype_list = line[8:]
			ref_hom_c = genotype_list.count("00")
			het_c = genotype_list.count("01")
			alt_hom_c = genotype_list.count("11")
			non_par_out.write(str(ref_hom_c)+"\t"+str(het_c)+"\t"+str(alt_hom_c)+"\n")
		elif line[0] == "superscaffold36":
			genotype_list = line[8:]
			ref_hom_c = genotype_list.count("00")
			het_c = genotype_list.count("01")
			alt_hom_c = genotype_list.count("11")
			scaf36_out.write(str(ref_hom_c)+"\t"+str(het_c)+"\t"+str(alt_hom_c)+"\n")
		elif line[0] == "superscaffold62":
			genotype_list = line[8:]
			ref_hom_c = genotype_list.count("00")
			het_c = genotype_list.count("01")
			alt_hom_c = genotype_list.count("11")
			scaf62_out.write(str(ref_hom_c)+"\t"+str(het_c)+"\t"+str(alt_hom_c)+"\n")
		else: # If the scaffold is in the non-recombining section of Z in females:
			genotype_list = line[8:]
			ref_hom_c = genotype_list.count("00")
			het_c = genotype_list.count("01")
			alt_hom_c = genotype_list.count("11")
			par_out.write(str(ref_hom_c)+"\t"+str(het_c)+"\t"+str(alt_hom_c)+"\n")

