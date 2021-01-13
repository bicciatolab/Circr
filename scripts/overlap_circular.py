#!/share/apps/bin/python
#run intersectBed like this:
#intersectBed -a /projects/seq-work/user/all/annotations/mouse/EnsemblGene-67.TR.gtf -b DP_2_circ.bed -wo > outputname.bed

#python overlap_circular.py bedfile gtffile limitlength
#bedfile ... output of bedtools intersect
#gtffile ... annotation file from ensembl
#limitlength ... discard circulars longer than x (use really high value to have all in the output)
#e.g.
#python overlap_circular.py bedfile.bed  /projects/seq-work/user/all/annotations/mouse/EnsemblGene-67.TR.gtf 500


########## SCRIPT MODIFIED by JIMMY   27020220

from collections import defaultdict
from functools import partial
from operator import itemgetter
from re import compile
from sys import argv
#from types import IntType

class GTF(object):
	def __init__(self, filename):
		self.__filename = filename
		self.__gtfdict = defaultdict(partial(defaultdict, list))
		self.__genreg = compile('ENSMUSG.[0-9]+') # MODIFIED ORIGINAL: compile('ENSG.[0-9]+')
		self.__transreg = compile('ENSMUST.[0-9]+') # MODIFIED ORIGINAL: compile('ENST.[0-9]+')
		self.__exonreg = compile('exon_number "([0-9]+)"') # MODIFIED ORIGINAL compile('exon_number "[0-9]+"')
		self.__genenamereg = compile('gene_name "(.*?)"')
		self.__transnamereg = compile('transcript_name "(.*?)"')
	def read_gtf(self):
		with open(self.__filename) as filein:
			for line in filein:
				linus = line.rstrip('\n').split('\t')
				if linus[2] == 'exon':
					last = linus[8]
					#print(last) # ADDED
					chrom, start, end, strand = linus[0], int(linus[3]) - 1, int(linus[4]), linus[6]
					geneid = self.__genreg.search(last).group()
					transid = self.__transreg.search(last).group()
					exonnr =  self.__exonreg.search(last).group(1)
					genename = self.__genenamereg.search(last).group(1)
					transname = self.__transnamereg.search(last).group(1)
					self.__gtfdict[(geneid, transid)][exonnr] = [(chrom, strand), (start, 's'), (end, 'e')]
	def get_gtfdict(self):
		return self.__gtfdict
	gtfdict = property(get_gtfdict)

class Overlap(object):
	def __init__(self, overlapfile, gtfdict, limitlength):
		self.__overlapfile = overlapfile
		self.__genreg = compile('ENSMUSG.[0-9]+') # ORIGINAL compile('ENSG.[0-9]+')
		self.__transreg = compile('ENSMUST.[0-9]+') # ORIGINAL compile('ENST.[0-9]+')
		self.__exonreg = compile('exon_number "([0-9]+)"') # MODIFIED ORIGINAL compile('exon_number "[0-9]+"')
		self.__genenamereg = compile('gene_name "(.*?)"')
		self.__transnamereg = compile('transcript_name "(.*?)"')
		self.__dictus = defaultdict(partial(defaultdict, set))
		self.__gtfdict = gtfdict
		self.__overlaplist = []
		self.__limitlength = limitlength

	def parse_overlap(self):
		with open(self.__overlapfile) as filein:
			for line in filein:
				linus = line.rstrip('\n').split('\t')
				geneinfo = linus[8]
				#only take lines from exon coords
				if linus[2] != 'exon':
					continue
				#get information for the gene/transcript
				#ensembl gtf -> start - 1
				gchrom, gstart, gend, gstrand = linus[0], int(linus[3]) - 1, int(linus[4]), linus[6]
				geneid = self.__genreg.search(geneinfo).group()
				transid = self.__transreg.search(geneinfo).group()
				exonnr =  int(self.__exonreg.search(geneinfo).group(1))
				genename = self.__genenamereg.search(geneinfo).group(1)
				transname = self.__transnamereg.search(geneinfo).group(1)
				#get information for the circular
				#bed format
				cchrom, cstart, cend, cname, cstrand = linus[9], int(linus[10]), int(linus[11]), linus[12], linus[14]
				#how much do they overlap
				overlap = int(linus[-1])
				if abs(cend - cstart) >= limitlength:
					continue
				self.__dictus[(geneid, transid, cname)]['chrom'].add((gchrom, cchrom))
				self.__dictus[(geneid, transid, cname)]['strand'].add((gstrand, cstrand))
				self.__dictus[(geneid, transid, cname)]['ccoord'].add((cstart, cend))
				self.__dictus[(geneid, transid, cname)]['name'].add((genename, transname))
				self.__dictus[(geneid, transid, cname)][exonnr] = (gstart, gend, overlap)

	def get_information_from_dictus(self):
		self.__overlaplist = ['\t'.join(['chromosom', 'geneid', 'genename', 'transcriptid', 'transcriptname', 'firstexon_start', 'lastexon_end', 'totalexon', 'overlapexon', 'exonlist', 'totallength', 'overlaplength', 'genestrand', 'circularid', 'circularstrand', 'circularstart', 'circularend', 'overlaptyp']) + '\n']
		for geneid, transid, cname in sorted(self.__dictus, key = itemgetter(2, 0, 1)):
			keytuple = (geneid, transid, cname)
			chrom = self.__dictus[keytuple]['chrom'].pop()[0]
			ccord = self.__dictus[keytuple]['ccoord'].pop()
			cstart, cend = ccord[0], ccord[1]
			strand = self.__dictus[keytuple]['strand'].pop()
			gstrand, cstrand = strand[0], strand[1]
			names = self.__dictus[keytuple]['name'].pop()
			genename, transcriptname = names[0], names[1]
			exonlist = sorted([i for i in self.__dictus[keytuple] if isinstance(i, int)]) ## prima 'int' era 'IntType'
			length = len(exonlist)

			if gstrand == '-':
				exonlist = exonlist[::-1]

			overlap_sum = sum([self.__dictus[keytuple][i][2] for i in exonlist])
			totalexon, totallength = self.get_total_exon_length(geneid, transid)

			firstexon_start, lastexon_end = self.__dictus[keytuple][exonlist[0]][0], self.__dictus[keytuple][exonlist[-1]][1]
			exonliststring = ','.join([str(i) for i in sorted(exonlist)])
			scstart, scend = cstart, cend

			if cstart == firstexon_start and cend == lastexon_end:
				if length == totalexon:
                                    overlapid = 'TTO'
                                    self.__overlaplist.append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\n'.format(chrom, geneid, genename, transid, transcriptname, firstexon_start, lastexon_end, totalexon, length, exonliststring, totallength, overlap_sum, gstrand, cname, cstrand, scstart, scend, overlapid))
				else:
                                    overlapid = 'PTO'
                                    self.__overlaplist.append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\n'.format(chrom, geneid, genename, transid, transcriptname, firstexon_start, lastexon_end, totalexon, length, exonliststring, totallength, overlap_sum, gstrand, cname, cstrand, scstart, scend, overlapid))
			

	def get_total_exon_length(self, gene, trans):
		exon = len(self.__gtfdict[(gene, trans)])
		length = 0
		for exonnr, entry in self.__gtfdict[(gene, trans)].items(): ### prima era iteritems
			length += abs(entry[2][0] - entry[1][0]) + 1
		return exon, length

	def circular_overlap_genes(self):
		self.__gene_circlist = ['\t'.join(('geneid', 'transcriptid', 'transcriptstrand', 'circularid', 'circularstrand', 'chromosom', 'circularstart', 'circularend')) + '\n']
		circset = set([i[2] for i in self.__dictus])
		for circular in circset:
			#get all genes/transcript which overlap with the circular
			genlist = [i for i in self.__dictus if i[2] == circular]
			#check if it overlaps more than one gene
			uniquegenes = set([i[0] for i in genlist])
			if len(uniquegenes) != 1:
				for entry in sorted(genlist, key = itemgetter(0, 1)):
					strand = list(self.__dictus[entry]['strand'])[0]
					chrom = list(self.__dictus[entry]['chrom'])[0]
					ccord = list(self.__dictus[entry]['ccoord'])[0]
					self.__gene_circlist.append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(entry[0], entry[1], strand[0], entry[2], strand[1], chrom[0], ccord[0], ccord[1]))

	def write_back(self, where, what):
		with open(where, 'w') as fileout:
			fileout.writelines(what)

	def get_overlaplist(self):
		return self.__overlaplist

	def get_genecirclist(self):
		return self.__gene_circlist

	overlaplist = property(get_overlaplist)
	genecirclist = property(get_genecirclist)

if __name__ == '__main__':
	bedfile = argv[1]
	gtffile = argv[2]
	limitlength = int(argv[3])
	ginst = GTF(gtffile)
	ginst.read_gtf()
	o = Overlap(bedfile, ginst.gtfdict, limitlength)
	o.parse_overlap()
	o.circular_overlap_genes()
	o.get_information_from_dictus()
	o.write_back(bedfile[:-3]+'overlap.txt', o.overlaplist)
	o.write_back(bedfile[:-3]+'gene_circ.txt', o.genecirclist)
