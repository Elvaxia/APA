import sys
from collections import Counter, OrderedDict
from pyfaidx import Fasta, wrap_sequence
fa, pos, outfa, outdisstat = sys.argv[1:]
fasta = Fasta(fa)
def makecord(site):
	"""
	return up- and down-stream region based on site
	"""
	gene, chr, strand, utrst, site, utred = site.split(':')
	return (chr, int(site) - 51, int(site) +50, strand)
	# return (chr, int(site) - 51, int(site), strand)
with open(pos) as P, open(outfa,'w') as OUT, open(outdisstat, 'w') as DIS:
	tmp = {}
	basedis = {}
	DIS.write('\t'.join(['location','A','C','G','T']) + '\n')
	for i in P.readlines():
		if i.startswith('PA'):continue
		#this is use the python spice method, so must substract one
		id, st, ed, strand = makecord(i.strip())
		seqname = id+':'+str(st)+'-'+str(ed)+':'+strand
		sequence = fasta[id][st:ed]
		line_len = fasta.faidx.index[id].lenc
		if strand == '-':
			sequence = sequence.complement
			sequence = sequence.reverse
			OUT.write('>' + seqname + '\n')
			for line in wrap_sequence(line_len, sequence.seq):
				OUT.write(line)
		else:
			OUT.write('>' + seqname + '\n')
			for line in wrap_sequence(line_len, sequence.seq):
				OUT.write(line)
		for index, base in enumerate(list(sequence.seq)):
			if index not in basedis:
				basedis[index] = [base]
			else:
				basedis[index].append(base)

	for index, base in basedis.items():
		c = Counter(base)
		total = sum(c.values())
		percent = {key: value/total for key, value in c.items()}
		try:
			del percent['N']
		except KeyError:
			pass
		percent = OrderedDict(sorted(percent.items()))
		DIS.write('\t'.join([str(index), '\t'.join(list(map(str, percent.values())))]) + '\n')
