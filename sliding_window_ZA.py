#!/usr/bin/python

import sys


"""
Written by Homa Papoli - August 2018
It takes a file containing scaffold name and their length,
with window size and step and outputs non-overlapping windows
of the given size. The output is in bed format so:
scaffold1 0 3
scaffold1 3 8

The first interval for scaffold1 contains bases 0,1,2
The second interval for scaffold1 contains bases 3,4,5,6,7

Input: length_f.txt
File containing two fields, scaffold name and scaffold length, for example:
NW_013188825.1 555

winSize is the window size which is an integer 
step is the step size which is an integer

The function ignores scaffolds shorter than the window size. 

Run as below:
python sliding_window_ZA.py length_f.txt winSize[int] step[int] > windows.txt
Example:
python sliding_window_ZA.py Anser_cygnoides.scaffold.length 100000 100000 > Anser_cygnoides.windows
"""


seq_l = open(sys.argv[1], 'r')
winSize = int(sys.argv[2])
step = winSize


# Function for sliding window
def slidingWindow(sequence_l, winSize, step):
	""" Returns a generator that will iterate through
	the defined chunks of input sequence. Input
	sequence must be iterable."""

	# Verify the inputs
	if not ((type(winSize) == type(0)) and (type(step) == type(0))):
		raise Exception("**ERROR** type(winSize) and type(step) must be int.")
	if step > winSize:
		raise Exception("**ERROR** step must not be larger than winSize.")
	if winSize > sequence_l:
		pass
	#winSize = sequence_l
	#numOfChunks = ((int(sequence_l-winSize)/step))+1
	#numOfChunks = int(numOfChunks)
	#else:
	# Pre-compute number of chunks to emit
	numOfChunks = ((int(sequence_l-winSize)/step))+1
	numOfChunks = int(numOfChunks) 

	for i in range(0, numOfChunks*step, step):
		yield i,i+winSize


for line in seq_l:
	sequence_l = int(line.strip("\n").split("\t")[1])
	if winSize > sequence_l:
		print(line.strip("\n").split("\t")[0]+"\t"+"0"+"\t"+str(sequence_l))
	else:
		for index, item in enumerate(slidingWindow(sequence_l, winSize, step)):
			print(line.strip("\n").split("\t")[0]+"\t"+str(item[0])+"\t"+str(item[1]))
			# For the last window when sequence_l%winSize != 0
			if (sequence_l-item[1]) < winSize and (sequence_l%winSize) != 0:
				print(line.strip("\n").split("\t")[0]+"\t"+str(item[1])+"\t"+str(sequence_l))

