#!/usr/bin/env python
'''
converts transcript gtf into exon read_gtf
be aware single exon transcript are lost!
input
argv[1] = gtf file
argv[2] = regex for transcript e.g. ENSMUST[0-9]+
'''


from re import compile
from re import search

from sys import argv

from collections import defaultdict
from functools import partial

def read_gtf(filename, regex):
	transdict = defaultdict(partial(defaultdict, int))
	stranddict = {}
	with open(filename) as filein:
		for line in filein:
			linus = line.rstrip('\n').split('\t')
			typus = linus[2]
			if typus != 'exon':
				continue
		
			transid = regex.search(linus[8]).group(0)
			es, ee = int(linus[3]), int(linus[4])
			exnr = int([i for i in linus[8].split(';') if 'exon_number'  in i ][0].split('"')[1])
			transdict[transid][exnr] = sorted((es, ee))
			stranddict[transid] = (linus[0], linus[6])
	return transdict, stranddict

def readinput(filename, transdict, stranddict):
	with open(filename) as filein:
		for line in filein:
			linus = line.split('\t')
			trans = linus[1]
			circ=linus[3]
			exon = int(linus[2])
			coords = transdict[trans][exon]
			strchr = stranddict[trans]
			print('{0}\t{1}\t{2}\t{4}\t{3}'.format(strchr[0], coords[0], coords[1], strchr[1], circ))


if __name__ == '__main__':
	regex = compile(argv[2])
	transdict, stranddict= read_gtf(argv[1], regex)
	inputfile = argv[3]
	readinput(inputfile, transdict, stranddict)
	