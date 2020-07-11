#!/usr/bin/env python
# encoding: utf-8
import sys
from itertools import islice
from pyfaidx import Fasta
import logging


logging.basicConfig(level=logging.INFO,format="%(levelname)s:%(asctime)s:%(message)s")

def binarySearch(structure_list, gene_region):
	"""
	Binary to restructure the motif CDS region

	binarySearch([1233,2233], 1000): (-1,0,False)
	binarySearch([1233,2233], 1300): (0,1,False)
	binarySearch([1233,2233], 1233): (0,0,True)
	binarySearch([1233,2233], 2400): (1,2,False)
	"""
	gene_region = int(gene_region)
	first = 0
	last = len(structure_list)-1
	found = False

	while first<=last and not found:
		pos = 0
		mid = (first + last)//2
		if structure_list[mid] == gene_region:
			pos = mid
			found = True
		else:
			if gene_region < int(structure_list[mid]):
				last = mid-1
			else:
				first = mid+1
	if found:
		return (pos, pos, found)
	else:
		return (last, first, found)

def GenerateDB(DBfile):
	"""
	recieve the polyA DB database, and return a DIC which contains the PA sites split by chromosome names.
	"""
	DBDIC = {}
	with open(DBfile) as IN:
		for i in IN.readlines():
			i = i.strip().split('\t')
			chrom, site = i[:2]
			if chrom not in DBDIC:
				DBDIC[chrom] = [site]
			else:
				DBDIC[chrom].append(site)
	return DBDIC
def window(seq, n=2):
	"""
    Return a sliding window over a list
    """
	it = iter(seq)
	result = list(islice(it, n))

	if len(result) == n:
		yield result
	for elem in it:
		result = result[1:] + list((elem,))
		yield result

def CalRatio(seq, strand):
	"""Receive a seq string, and return the biggest ratio.
	"""
	# if strand == '-':
	# 	A_ = seq.count("A")
	# else:
	# 	A_ = seq.count("T")
	# try:
	# 	Allcount = len(seq)
	# 	return A_/Allcount
	# except:
	# 	#exclude the "N"
	# 	return 0
	A = "TTTTTT" if strand == '+' else "AAAAAA"
	if A in seq:
		return False
	else:
		return True

def CalAorT(line, fasta):
	# info = line.strip().split('\t')
	chrom, st, ed = line[:3]
	strand = "+" if line[-1] == "-" else "-"
	peak = int((int(st) + int(ed)) / 2 )
	#这个地方发现有些上下游刚好的polyA的中间，所以这里计算比例是否得考虑更换
	poly = "TTTTTT" if strand == '+' else "AAAAAA"
	# if strand == '-':
	# 	seq =  fasta[chrom][peak - 15: peak + 15]
	# else:
	seq = fasta[chrom][peak - 40:peak + 40]

	if poly in str(seq):
		return False
	else:
		return True
	
def main():
	fa, bedfile, DBfile, fileout = sys.argv[1:]
	fasta = Fasta(fa)
	PROCOUNT  = 0
	INFOMAX   = 100000
	logging.info('Process the polyDB information')
	PAdb	  = GenerateDB(DBfile)

	logging.info("Start process region file")

	with open(bedfile) as bed, open(fileout, 'w') as O:
		for i in bed.readlines():
			PROCOUNT += 1
			info = i.strip().split('\t')
			##check the intersect with known database


			if PROCOUNT % INFOMAX == 0:
				O.flush()
				logging.info("Parsed {} regeions.".format(PROCOUNT))
			intersect = binarySearch(PAdb[info[0]],info[1])
			if intersect[-1]:
				O.write(i)
			else:
				try:
					if info[2] > PAdb[info[0]][intersect[1]]:
						O.write(i)
					else:
						if CalAorT(info, fasta):
							O.write(i)
						else:
							# print(i.strip())
							pass
				except IndexError:
					if CalAorT(info, fasta):
						O.write(i)
					# else:
						# print(i.strip())
						# pass
				

if __name__ == '__main__':
	main()
