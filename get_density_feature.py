#!/usr/bin/python
from __future__ import division
import sys
import pprint
from fasta import readfasta

#************************************************
# Written by Hom Papoli - 18 Nov 2018
# Calculates window and feature base and SNP count
# ./get_density_het_window.py overlap.bed seq.fasta
#************************************************
# Define GC_count_f function
def GC_count_f(sequence):
	"""Return heterozygosity from
		any sequence as a string """
	# The first thing is to test the input
	if not isinstance(sequence, str):
		raise Exception("Sequence is not a string")
	homozygote='ATCG'
	GCs='GC'
	homozygote_count = len([base.upper() for base in sequence if base.upper() in homozygote])
	GC_count = len([base.upper() for base in sequence if base.upper() in GCs])
	return homozygote_count, GC_count
#************************************************
# Read bed and fasta input files
#************************************************
overlapping = open(sys.argv[1], "r")
fasta_seq = open(sys.argv[2], "r")
out = open(sys.argv[3], "w")
#************************************************
# Read fasta file in dictionary
fasta_dict = readfasta(fasta_seq)
#************************************************
#************************************************
# Read overlap file in dictionary
#************************************************
overlapping_list = [] # Create a list to keep order when printing
overlapping_dict = {}
multi_window_intervals = {}
for line in overlapping:
	line = line.strip("\n").split("\t")
	key = line[0]+":"+line[1]+":"+line[2]
	value = [int(line[4]), int(line[5]), int(line[6])]
	if key in overlapping_dict.keys():
		overlapping_dict[key].append(value)
	else:
		overlapping_dict[key] = [value]
	if not key in overlapping_list:
		overlapping_list.append(key)

#pp = pprint.PrettyPrinter(indent=4)
#pp.pprint(overlapping_dict)

#************************************************
# Calculate base count and feature
# density for every window.
#************************************************
header = ["Scaffold", "Window_start", "Window_end", "Window_Base_count", "Window_GC_count", "Feat_start", "Feat_end", "Feat_Base_count", "Feat_GC_count"]
#print("\t".join(header))
out.write("\t".join(header)+"\n")
for key in overlapping_list:
	print(key)
	bases = []
	gcs = []
	scaffold = key.split(":")[0]
	winstart = int(key.split(":")[1])
	winend = int(key.split(":")[2])
	for interval in overlapping_dict[key]:
		featstart = interval[0]
		featend = interval[1]
		sequence = fasta_dict[scaffold][featstart:featend]
		Base_count = list(GC_count_f(sequence))[0]
		GC_count = list(GC_count_f(sequence))[1]
		bases.append(Base_count)
		gcs.append(GC_count)
	sequence_widnow = fasta_dict[scaffold][winstart:winend]
	Base_count_window = list(GC_count_f(sequence_widnow))[0]
	GC_count_window = list(GC_count_f(sequence_widnow))[1]
	info = [scaffold, winstart, winend, Base_count_window, GC_count_window, featstart, featend, sum(bases), sum(gcs)]
	out.write("\t".join(str(i) for i in info)+"\n")

overlapping.close()
fasta_seq.close()
