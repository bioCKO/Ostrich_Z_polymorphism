#!/usr/bin/python
from __future__ import division
import sys
import operator as op
import math
#from functools import reduce

# **********************************
# Written by Homa Papoli - Nov 2018
# This script gets two population
# data and outputs pi and dxy 
# **********************************


# **********************************
# Function denfinitions
# **********************************
# Function for N choose r
# **********************************
#def ncr(n, r):
#    r = min(r, n-r)
#    numer = reduce(op.mul, xrange(n, n-r, -1), 1)
#    denom = reduce(op.mul, xrange(1, r+1), 1)
#    return numer//denom
def ncr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)
# **********************************
# Function for Pairwise differences
# **********************************
def pi_fun(nchr, allele1, allele2):
	P_diff = allele1*allele2
	NCR = ncr(nchr, 2)
	pi = round(P_diff/NCR, 4)
	return (P_diff, NCR, pi) 
# **********************************
# Function for dxy
# **********************************
def dxy_fun(nchr_pop1, nchr_pop2, allele1_pop1, allele1_pop2):
	p1 = allele1_pop1/nchr_pop1
	p2 = allele1_pop2/nchr_pop2
	dxy = round(p1*(1-p2) + p2*(1-p1), 4)
	return dxy
# **********************************
# Function for FST
# **********************************

# **********************************
# Open input
# **********************************
pop1 = open(sys.argv[1], "r")
# **********************************
# Open output for counts of number
# of pairwise differences and number
# of pairwise comparisons
# **********************************
counts_out = open(sys.argv[2], "w")
#pop2 = open(sys.argv[2], "r")
#pop1_pi = open(sys.argv[3], "w")
#pop2_pi = open(sys.argv[4], "w")
#dxy_out = open(sys.argv[5], "w")
# **********************************
# Read pop1 and pop2 into dictionary
# **********************************
#pop1_dict = {}
#print("Scaffold"+"\t"+"SNP_pos"+"\t"+"Pi_pop1"+"\t"+"Pi_pop2"+"\t"+"dxy_pop12")
for line in pop1:
	#print(line)
	line = line.strip("\n").split()
	#print(line[3], line[4].split(":")[0], line[5])
	if not line[0] == "scaffold":
	#	print(int(line[3]), line[4])#int(line[4].split(":")))#[1]), int(line[5].split(":")[1]))
		pop1_snps_counts = list(pi_fun(int(line[3]), int(line[4].split(":")[1]), int(line[5].split(":")[1])))[0]
		pop1_base_counts = list(pi_fun(int(line[3]), int(line[4].split(":")[1]), int(line[5].split(":")[1])))[1]
		pop1_pi = list(pi_fun(int(line[3]), int(line[4].split(":")[1]), int(line[5].split(":")[1])))[2]
		#print(pop1_pi)
		#pop2_pi = pi_fun(int(line[9]), int(line[10]), int(line[11]))
		#dxy_pop12 = dxy_fun(int(line[3]), int(line[9]), int(line[4]), int(line[10]))
		print(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+str(pop1_pi)+"\t"+"\t".join(line[6:]))#+"\t"+line[7]+"\t"+line[8]+"\t"+line[9]+"\t"+line[10]+"\t"+line[11]+"\t"+line[12]+"\t"+line[13]+"\t"+line[14]+"\t"+line[15]+"\t"+line[16]+"\t"+line[17])#+"\t"+str(pop2_pi)+"\t"+str(dxy_pop12))
		counts_out.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+str(pop1_pi)+"\t"+str(pop1_snps_counts)+"\t"+str(pop1_base_counts)+"\t"+"\t".join(line[6:])+"\n")
	#else:
	#	print("\t".join(line[0:6])+"\t"+"pi"+"\t"+"\t".join(line[6:]))


# 		print(line[0], line[1], line[2], line[3], line[4].split(":")[1], line[5].split(":")[1])
# 		#key, value = line[0]+"_"+line[1], line[3:]
# 		#if key in pop1_dict.keys():
# 		#	raise Exception
# 		#else:
# 		#	pop1_dict[key] = value

# pop2_dict = {}
# for line in pop2:
# 	line = line.strip("\n").split("\t")
# 	if not line[0] == "CHROM":
# 		key, value = line[0]+"_"+line[1], line[3:]
# 		if key in pop2_dict.keys():
# 			raise Exception	
# 		else:
# 			pop2_dict[key] = value
# #print(pop1_dict, pop2_dict)
# # **********************************
# # Calculate pi for each SNP
# # **********************************
# # pop1
# # **********************************
# for snp in pop1_dict.keys():
# 	nchr = int(pop1_dict[snp][0])
# 	allele1 = int(pop1_dict[snp][1].split(":")[1])
# 	allele2 = int(pop1_dict[snp][2].split(":")[1])
# 	pop1_pi.write(snp.split("_")[0]+"\t"+snp.split("_")[1]+"\t"+str(pi_fun(nchr, allele1, allele2))+"\n")
# # **********************************
# # pop2
# # **********************************
# for snp in pop2_dict.keys():
# 	nchr = int(pop2_dict[snp][0])
# 	allele1 = int(pop2_dict[snp][1].split(":")[1])
# 	allele2 = int(pop2_dict[snp][2].split(":")[1])
# 	pop2_pi.write(snp.split("_")[0]+"\t"+snp.split("_")[1]+"\t"+str(pi_fun(nchr, allele1, allele2))+"\n")
# # **********************************
# # Calculate dxy for each SNP
# # **********************************
# for snp in pop1_dict.keys():
# 	nchr1 = int(pop1_dict[snp][0])
# 	nchr2 = int(pop2_dict[snp][0])
# 	allele1_1 = int(pop1_dict[snp][1].split(":")[1])
# 	allele1_2 = int(pop2_dict[snp][1].split(":")[1])
# 	dxy_out.write(snp.split("_")[0]+"\t"+snp.split("_")[1]+"\t"+str(dxy_fun(nchr1, nchr2, allele1_1, allele1_2))+"\n")
