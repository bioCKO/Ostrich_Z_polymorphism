#!/usr/bin/python
from __future__ import division
import sys
import numpy as np
from operator import itemgetter
from itertools import chain
from datetime import datetime
from scipy import stats
import math

"""
Written by Homa Papli - Jan 2019
It outputs confidence interval of the mean for a given resampling distribution. 
"""
f1 = open(sys.argv[1], "r")

f1_dict = {}
f1_list = []
for line in f1:
	line = line.strip("\n").split()
	key = line[0]+"_"+line[12]+"_"+line[13]
	value = line
	if not key in f1_list:
		f1_list.append(key)
	if key in f1_dict.keys():
		f1_dict[key].append(value)
	else:
		f1_dict[key] = [value]

resample_dict = {}
for key in f1_list:
	pi_count = []
	for element in f1_dict[key]:
		#print(element, element[6])
		pi_count.append(float(element[6]))
	pi_count_key = sum(pi_count)
	if key.split("_")[0] == "superscaffold54":
		if int(key.split("_")[1]) > 16379243: # Removing the non-Z linked part of superscaffold54
			pass
		else:
			if key in resample_dict.keys():
				resample_dict[key].append([pi_count_key, int(element[26])])
			else:
				resample_dict[key] = [[pi_count_key, int(element[26])]]
	else:
		if key in resample_dict.keys():
			resample_dict[key].append([pi_count_key, int(element[26])])
		else:
			resample_dict[key] = [[pi_count_key, int(element[26])]]

#print(resample_dict)

for i in range(0, 99):	
	#print(i)	
	update_dict = {}
	for key in f1_dict.keys():
		if key.split("_")[0] == "superscaffold54":
			if int(key.split("_")[1]) > 16379243: # Removing the non-Z linked part of superscaffold54
				pass
			else:
				# resample the list of a given window of a given scaffold
				list_len = len(f1_dict[key])
				rn = list(np.random.randint(low=0, high=list_len, size=list_len))
				overlap = list(itemgetter(*rn)(f1_dict[key]))
				update_dict[key] = overlap
		else:
			list_len = len(f1_dict[key])
			rn = list(np.random.randint(low=0, high=list_len, size=list_len))
			overlap = list(itemgetter(*rn)(f1_dict[key]))
			update_dict[key] = overlap			

#for key in update_dict.keys():
#	if "superscaffold36" in key:
#		print(key, update_dict[key])

	for key in f1_list:
		#print(key)
		if key in update_dict.keys():
			pi_count = []
			for element in update_dict[key]:
				pi_count.append(float(element[6]))
				#print(pi_count)
			pi_count_key = sum(pi_count)
			#print(pi_count_key)
			if key.split("_")[0] == "superscaffold54":
				if int(key.split("_")[1]) > 16379243: # Removing the non-Z linked part of superscaffold54
					pass
				else:
					if key in resample_dict.keys():
						resample_dict[key].append([pi_count_key, int(element[26])])
					else:
						resample_dict[key] = [[pi_count_key, int(element[26])]]
			else:
				if key in resample_dict.keys():
					resample_dict[key].append([pi_count_key, int(element[26])])
				else:
					resample_dict[key] = [[pi_count_key, int(element[26])]]
#print(len(resample_dict['superscaffold69-1_1000000_2000000']))

for key in f1_list:
	if key in resample_dict.keys():
		#print(resample_dict[key])
		list_len = len(f1_dict[key])
		#print(list_len)
		df_N = list_len-1
		snps = []
		bases = []
		het = []
		for element in resample_dict[key]:
			snps.append(element[0])
			bases.append(element[1])
			het.append(element[0]/element[1])
		mean_of_resampling = np.nansum(snps)/np.nansum(bases)
		#print(mean_of_resampling)
		t_student = stats.t(df=df_N).ppf((0.025, 0.975))
		std = np.std(het)
		CI_lower = mean_of_resampling - list(t_student)[1]*(std)#/math.sqrt(intergenic_len))
		CI_upper = mean_of_resampling + list(t_student)[1]*(std)
		print(key.split("_")[0]+"\t"+key.split("_")[1]+"\t"+key.split("_")[2]+"\t"+str(mean_of_resampling)+"\t"+str(CI_lower)+"\t"+str(CI_upper)+"\t"+str(std))


		



