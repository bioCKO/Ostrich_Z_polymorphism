#!/usr/bin/python
import sys
import collections

#########################################################################################################################
# Written by Homa Papoli - 26 May 2016  
# Updated 16 Feb 2017   
# Update 30 Sept 2018                                                                             #
# This script takes a gff file as an input and outputs the coordinates                                                  #
# for intergenic regions, intron or exons based on the term determined                                                  #
# in the argument.                               																	    #
# For intergenic, run:                                                  									    #
# Usage: ./extract_from_gff.py gff_file [intergenic/intron/CDS] [first/all/second_to_last] scaffold_length.txt > output #
# For CDS and intronic, run:																						    #
# Usage: ./extract_from_gff.py gff_file [intergenic/intron/CDS] [first/all/second_to_last] > output 
# values for CDS_coordinates dictionary differ between different gtf/gff formats and should be adapted accordingly.
#########################################################################################################################

#*****************************************************************************
# Specify the inputs
#*****************************************************************************
gff = open(sys.argv[1], 'r')
element = sys.argv[2] # intergenic, intron or CDS
number = sys.argv[3] # if first intron/exon, enter 'first', if all introns/exons, enter 'all', 
					# if second to last introns/exons, enter 'second_to_last'
					# for intergenic, the third argument is always 'all'
					
# For the intergenic regions, the start and end of the scaffolds should also be
# taken into account. For this, we need the scaffold lengths.					
# Read scaffold lengths as a dictionary:
if element == "intergenic":
	scaf = open(sys.argv[4], "r")
	scaf_dict = {}
	for line in scaf:
		line = line.strip("\n").split("\t")
		scaffold, length = line[0], int(line[1])
		scaf_dict[scaffold] = length
    #print scaf_dict
#*****************************************************************************

#*****************************************************************************
# Read the mRNA and CDS coordinates into dictionaries
#*****************************************************************************
mRNA_coordinates = {}
CDS_coordinates = {}
region_coordinates = {} # Some scaffolds do not contain any gene which were ignored previously. 
for line in gff:
	if not line.startswith("#"):
		line = line.strip('\n').split('\t')
		if line[2] == 'mRNA': # In some formats, it's mRNA and in others, it's gene. 'mRNA': # if the third column of the gff file is mRNA
			#print(line)
			key, value = line[0], line[3:5] # set the scaffold name as key and the start:stop as value
			if key in mRNA_coordinates.keys(): # if the scaffold name is already present in the mRNA_coordinates dictionary
				mRNA_coordinates[key].append(value) # append the new coordinates to the list
			else:
				mRNA_coordinates[key] = [value] # otherwise set a new coordinate list to the new key
		elif line[2] == 'CDS': # if the third column of the gff file is CDS
            #if 'Genbank' in line[8].split(";")[2]:
                #print(line[8])
			key, value = line[0]+'_'+line[8].replace(';', '').replace('Parent=', ''), [line[3], line[4], line[6], line[7]] # set the key in the form scaffold_genename and the start:stop as value 
			if key in CDS_coordinates.keys(): # if the third column of the gff file is CDS
				CDS_coordinates[key].append(value) 
			else:
				CDS_coordinates[key] = [value]
		elif line[2] == 'region':
			key, value = line[0], line[3:5]
			if key in region_coordinates.keys():
				region_coordinates[key].append(value)
			else:
				region_coordinates[key] = [value]
# NW_013185661.1_cds3936_XP_013027331.1	3383542	3383711	-
# NW_013185654.1  Gnomon  CDS     6       296     .       +       0       ID=cds0;Parent=rna0;Dbxref=GeneID:106029314,Genbank:XP_013029950.1;Name=XP_013029950.1;gbkey=CDS;gene=LOC1060
no_gene_scaffold = {}
for scaf in region_coordinates.keys():
	if not scaf in mRNA_coordinates.keys():
		no_gene_scaffold[scaf] = region_coordinates[scaf]
#print(no_gene_scaffold.keys())
#*****************************************************************************
# Read the intergenic coordinates into a dictionary
#*****************************************************************************			
#print(mRNA_coordinates['NW_010286248.1'])
#print(mRNA_coordinates['NW_010290423.1'])	
#print mRNA_coordinates
if element == 'intergenic' and number == 'all': # element is the second argument of the command line, if it is intergenic
	intergenic = {} # define a dictionary called intergenic
	for scaffold in mRNA_coordinates.keys(): # look for the scaffold in the keys of mRNA_coordinates dictionary
		if len(mRNA_coordinates[scaffold]) == 1:
			#intergenic[scaffold] = mRNA_coordinates[scaffold]
			print(scaffold+'\t'+'0'+'\t'+str(mRNA_coordinates[scaffold][0][0]))
			print(scaffold+'\t'+str(int(mRNA_coordinates[scaffold][0][1]))+'\t'+str(scaf_dict[scaffold]))
		else:
			for i in range(len(mRNA_coordinates.get(scaffold))): # calculate the length of each value of the mRNA_coordinates dictionary and store its range in a variable called i 
				#print(i, len(mRNA_coordinates['NW_010286248.1']))
				if i == len(mRNA_coordinates.get(scaffold)) - 1: # if i is equal to the length - 1 because i starts from 0 so if you have l = [1, 2, 3, 4], [i for i in range(len(l))] would be [0, 1, 2, 3]
					break # so the last element would have an index which is equal to the length of the list minus 1, if so, stop the loop
				else: # otherwise
					if scaffold in intergenic.keys(): # if scaffold exists among the keys in the intergenic dictionary
						intergenic[scaffold].append([mRNA_coordinates.get(scaffold)[i][1], mRNA_coordinates.get(scaffold)[i+1][0]]) # append to the list the intergenic coordinates
					else:
						intergenic[scaffold] = [[mRNA_coordinates.get(scaffold)[i][1], mRNA_coordinates.get(scaffold)[i+1][0]]]
	#print(intergenic['NW_010290423.1'])
	#print(intergenic)				
	#**********************************************
	# Print the intergenic coordinate to the output
	#**********************************************
	# In order to obtain the intergenic sequences from the gff file, the sequences
	# between the two mRNA coordinates were obtained. However, the start and end
	# sequence of each scaffold is also potentially part of the intergenic sequence.
	# The code below takes into account these segments as well.			
	for key, value in intergenic.iteritems():
		#print key, value
		#Example when there are more than 2 coordinate sets:
		#scaffold3301 [['27574', '38089'], ['38253', '40600'], ['40764', '72871']]
		#scaffold3301    0       27573
		#scaffold3301    27573   38089
		#scaffold3301    38252   40600
		#scaffold3301    40763   72871
		#scaffold3301    72871   73321
		if len(value) > 2:
			for element in value:
				if value.index(element) == 0:
					#print (key+'\t'+"0"+'\t'+str(int(mRNA_coordinates['NW_010290423.1'][0][0])-1))
					print (key+'\t'+"0"+'\t'+str(int(mRNA_coordinates[key][0][0])))
					#print (key+'\t'+"0"+'\t'+str(int(element[0])-1))
					print (key+'\t'+str(int(element[0])-1)+'\t'+element[1])
				elif value.index(element) == len(value)-1:
					if key in scaf_dict.keys():
						print (key+'\t'+str(int(element[0])-1)+'\t'+element[1])
						print (key+'\t'+mRNA_coordinates[key][-1][1]+'\t'+str(scaf_dict[key]))
				else:
					print (key+'\t'+str(int(element[0])-1)+'\t'+element[1])
#		scaffold54855 [['6377', '15720']]
#		scaffold54855   0       6376
#		scaffold54855   6376    15720
#		scaffold54855   15720   17673	
		elif len(value) == 1:
			for element in value:
				print (key+'\t'+"0"+'\t'+str(int(mRNA_coordinates[key][0][0])))
				print (key+'\t'+str(int(element[0])-1)+'\t'+element[1])
				if key in scaf_dict.keys():
					print (key+'\t'+mRNA_coordinates[key][-1][1]+'\t'+str(scaf_dict[key]))
#		scaffold8907 [['56060', '68023'], ['85306', '102588']]
#		scaffold8907    0       56059
#		scaffold8907    56059   68023
#		scaffold8907    85305   102588
#		scaffold8907    102588  126423
		elif len(value) == 2:
			for element in value:
				if value.index(element) == 0:
					print (key+'\t'+"0"+'\t'+str(int(mRNA_coordinates[key][0][0])))
					#print (key+'\t'+"0"+'\t'+str(int(element[0])-1))
					print (key+'\t'+str(int(element[0])-1)+'\t'+element[1])
				elif value.index(element) == 1:
					if key in scaf_dict.keys():
						print (key+'\t'+str(int(element[0])-1)+'\t'+element[1])
						print (key+'\t'+mRNA_coordinates[key][-1][1]+'\t'+str(scaf_dict[key]))
	
	# Print scaffolds with no genes on them
#	for i in no_gene_scaffold.keys():
#		print (i+'\t'+str(int(no_gene_scaffold[i][0][0])-1)+'\t'+str(int(no_gene_scaffold[i][0][1])))
# #*****************************************************************************
# # Read the intronic coordinates into a dictionary
# #*****************************************************************************		
if element == 'intron':	# if element is intron
	#print(CDS_coordinates)
	intronic = {} # define a dictionary called intron
	for scaffold in CDS_coordinates.keys(): # look for the scaffold in the keys of CDS_coordinates dictionary
		if len(CDS_coordinates[scaffold]) != 1: # if length of the value list is not equal to 1, because if it is equal to 1, it means there is only one CDS and no intron
			if CDS_coordinates[scaffold][0][2] == "+":
				for position1, position2 in zip(CDS_coordinates[scaffold], CDS_coordinates[scaffold][1:]): # now use the zip function which pairs every two elements in a list
					#print("zip", zip(CDS_coordinates[scaffold], CDS_coordinates[scaffold][1:]))
					if scaffold in intronic.keys():
						intronic[scaffold].append([position1[1], position2[0]])
					else:	
						intronic[scaffold] = [[position1[1], position2[0]]]
			else:
				for position1, position2 in zip(CDS_coordinates[scaffold], CDS_coordinates[scaffold][1:]): # now use the zip function which pairs every two elements in a list
					#print("zip", zip(CDS_coordinates[scaffold], CDS_coordinates[scaffold][1:]))
					if scaffold in intronic.keys():
						intronic[scaffold].append([position2[1], position1[0]])
					else:	
						intronic[scaffold] = [[position2[1], position1[0]]]

	#*********************************************
	# Print all intronic coordinates to the output
	#*********************************************				
	if number == 'all': # in case the second argument is intron, the third argument can take 'all', 'first' or 'second_to_last'
		for key, value in intronic.iteritems():
			for element in value:
				print (key+'\t'+str(int(element[0]))+'\t'+str(int(element[1])-1))
				
	#***************************************************
	# Print the first intronic coordinates to the output
	#***************************************************
	if number == 'first': # if it is first, it only ouputs the coordinates of the first intron
		for key, value in intronic.iteritems():
			if len(value)>1:
				print (key+'\t'+str(int(value[0][0]))+'\t'+str(int(value[0][1]-1)))
				
	#************************************************************
	# Print from the second to the last coordinates to the output
	#************************************************************				
	if number == 'second_to_last': # if it is second_to_last it outputs the intronic coordinates from the second to the last one
		for key, value in intronic.iteritems():
			if len(value)>1:
				for index, item in enumerate(value):
					if index != 0:
						print (key+'\t'+str(int(item[0]))+'\t'+str(int(item[1]-1)))
			elif len(value) == 1:
				print (key+'\t'+str(int(value[0][0]))+'\t'+str(int(value[0][1]-1)))
						
#************************************************
# Print all the CDS coordinates into a dictionary
#************************************************							
if element == 'CDS' and number =='all': # the same thing applies to the CDS, if it is all, it prints the coordinates for all CDS
	for key, value in CDS_coordinates.iteritems():
		for element in value:
			print (key+'\t'+str(int(element[0])-1)+'\t'+element[1]+'\t'+element[2])
			
#*************************************************
# Print the first CDS coordinate into a dictionary
#*************************************************		
elif element == 'CDS' and number == 'first': # if it is first, it prints the coordinates for the first CDS only
	for key, value in CDS_coordinates.iteritems():
		if len(value)>1:
			print (key+'\t'+str(int(value[0][0])-1)+'\t'+value[0][1])
			
#**************************************************************
# Print the second to the last CDS coordinate into a dictionary
#**************************************************************
elif element == 'CDS' and number == 'second_to_last': # if it is second_to_last, it prints the CDS coordinates for the second one to the last
	for key, value in CDS_coordinates.iteritems():
		if len(value)>1:
			for index, item in enumerate(value):
				if index !=0:
					print (key+'\t'+str(int(item[0])-1)+'\t'+item[1])
		elif len(value) == 1:
			print (key+'\t'+str(int(value[0][0])-1)+'\t'+value[0][1])
	
			
gff.close()
if element == "intergenic":
	scaf.close()
