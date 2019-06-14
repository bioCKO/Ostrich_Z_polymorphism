#!/usr/bin/python
from __future__ import division
import sys
import pprint
#from fasta import readfasta
#from het import heterozygosity

#************************************************
# Written by Hom Papoli - 22 Oct 2018
# Correct and re-write the overlapping list of 
# intervals for cases where an interval corresponds 
# to several windows. For example:
# NW_013185928.1  200000  300000  NW_013185928.1  266985  310105  33015
# NW_013185928.1  300000  400000  NW_013185928.1  266985  310105  10105
# After correction:
# NW_013185928.1	200000	300000	NW_013185928.1	266985	300000	33015
# Input:
# scaffold    Start    End    Scaffold    Start    End    Length
# NW_013185928.1	0	100000	NW_013185928.1	0	2082	2082
# NW_013185928.1	0	100000	NW_013185928.1	10047	10260	213
# NW_013185928.1	0	100000	NW_013185928.1	64687	65241	554
# NW_013185928.1	0	100000	NW_013185928.1	73623	75089	1466
# Run:
# ./get_sum_overlapping.py overlap.bed
#************************************************

#************************************************
# Read input
#************************************************
overlapping = open(sys.argv[1], "r")
#************************************************
# Read fasta file in dictionary
#fasta_dict = readfasta(fasta_seq)
#************************************************

#************************************************
# Read overlap file in dictionary
#************************************************
overlapping_dict = {}
multi_window_intervals = {}
for line in overlapping:
	line = line.strip("\n").split("\t")
	 # When there is no overlap between A and B, start and end are reported with "-1".
	if line[4] == '-1':
		line[3] = line[0]
		line[4] = '0'
		line[5] = '0'
		line[6] = '0'
		key = line[3]+":"+line[4]+":"+line[5]
		value = [int(line[1]), int(line[2]), int(line[4]), int(line[5]), int(line[6])]
		if key in overlapping_dict.keys():
			overlapping_dict[key].append(value)
		else:
			overlapping_dict[key] = [value]
	else: # In all other cases:
		# Value is [winstart, winend, featurestart, featureend, overlaplength]
		key = line[3]+":"+line[4]+":"+line[5]
		value = [int(line[1]), int(line[2]), int(line[4]), int(line[5]), int(line[6])]
		# Append such cases into the overlapping_dict
		if key in overlapping_dict.keys():
			overlapping_dict[key].append(value)
		else:
			overlapping_dict[key] = [value]
	
#print(overlapping_dict)	
#pp = pprint.PrettyPrinter(indent=4)
#pp.pprint(overlapping_dict)
#print(len(overlapping_list))
#print(len(overlapping_dict.keys()))

#************************************************
# Update the overlap dictionary
#************************************************
# Correct and re-write the overlapping list of 
# intervals for cases where an interval corresponds 
# to several windows. For example:
# NW_013185928.1  200000  300000  NW_013185928.1  266985  310105  33015
# NW_013185928.1  300000  400000  NW_013185928.1  266985  310105  10105
# After correction:
# NW_013185928.1	200000	300000	NW_013185928.1	266985	300000	33015
# In this case, new coordinates 
# for the features should be generated. 
#************************************************
# Define a new dictionary
update_dict = {} 
for key in overlapping_dict.keys():
	# in that case, it means there are several windows for one feature interval.
	if len(overlapping_dict[key]) > 1: 
		# this is the length of features put into one list
		new_value = [i[4] for i in overlapping_dict[key]]
		for index, item in enumerate(new_value): # loop through the feature length list
			if index == 0: # for the first length
				start = int(key.split(":")[1]) # start is the start of the feature
				end = int(key.split(":")[1]) + new_value[index] # end is the start of the featue plus the length of the feature
				# update_key is scaffold:winstart:winend
				update_key = key.split(":")[0]+":"+str(overlapping_dict[key][index][0])+":"+str(overlapping_dict[key][index][1])
				# update_value is [featurestart, featureend, featurelength]
				update_value = [start, end, end-start]
				if update_key in update_dict.keys():
					update_dict[update_key].append(update_value)
				else:
					update_dict[update_key] = [update_value]
			else: # for indexes larger than 0
				# start is the start of the feature plus the length of the previous, 
				# as start is then saved, this works.
				start = int(key.split(":")[1]) + sum([new_value[i] for i in range(0, index)])#new_value[index-1] 
				end = int(key.split(":")[1]) + sum([new_value[i] for i in range(0, index)]) + new_value[index]
				update_key = key.split(":")[0]+":"+str(overlapping_dict[key][index][0])+":"+str(overlapping_dict[key][index][1])
				update_value = [start, end, end-start]
				if update_key in update_dict.keys():
					update_dict[update_key].append(update_value)
				else:
					update_dict[update_key] = [update_value]		
	else: # If there is only one window for a given feature interval.
		update_key = key.split(":")[0]+":"+str(overlapping_dict[key][0][0])+":"+str(overlapping_dict[key][0][1])
		update_value = [overlapping_dict[key][0][2], overlapping_dict[key][0][3], overlapping_dict[key][0][4]]
		if update_key in update_dict.keys():
			update_dict[update_key].append(update_value)
		else:
			update_dict[update_key] = [update_value]

#pp = pprint.PrettyPrinter(indent=4)
#pp.pprint(update_dict)

for key in update_dict.keys():
	for item in update_dict[key]:
		l = [key.split(":")[0], key.split(":")[1], key.split(":")[2], key.split(":")[0], item[0], item[1], item[2]]
		print("\t".join(str(i) for i in l))

overlapping.close()

