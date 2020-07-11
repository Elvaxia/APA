import sys
import os
import time
import pysam
import tables
import logging

from heapq import nlargest
from operator import itemgetter
from itertools import islice, groupby, chain, combinations, product
from collections import Counter, defaultdict, OrderedDict
from math import ceil

from operator import itemgetter
from multiprocessing import Pool

logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(asctime)s:%(message)s")


class ReadParser(object):
    def __init__(self, readline):
        """
        继承了pysam包的AlignedSegment这一类，直接搞。加入了想要的pA位点的判定
        :param readline:
        """
        self.readline = readline
        self.cigarinfo = self.readline.cigartuples
        self.sequence = self.readline.seq

    @property
    def softclip(self):
        """
        :return: soft clip sequence
        """
        if not self.readline.is_reverse:
            if self.cigarinfo[-1][0] != 4:
                return False
            else:
                if self._lastcigar == 1:
                    return False
                else:
                    tmpseq = self.sequence[-self._lastcigar:]
                    if not self.__checkpolyA(tmpseq):
                        return False
                    return True
        else:
            if self.cigarinfo[0][0] != 4:
                return False
            else:
                if self._lastcigar == 1:
                    return False
                else:
                    tmpseq = self.__rev_comseq(self.sequence[:self._lastcigar])
                    if not self.__checkpolyA(tmpseq):
                        return False
                    return True

    @property
    def _lastcigar(self):
        """
        :return: return cigar length of the last 3' direction
        """
        if not self.readline.is_reverse:
            return int(self.cigarinfo[-1][1])
        else:
            return int(self.cigarinfo[0][1])

    def __rev_comseq(self, seq):
        """
        :return: reverse compliment the sequence
        """
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
        return reverse_complement

    def __checkpolyA(self, seq):
        """
        Checking the soft clip reads, continuous A or not
        :return: True of False
        """

        def headseq(seq):
            """Return the top six sequence
            >>> Headseq('ATCTAAA')
            'ATCTAA'
            >>> Headseq('ATCT')
            'ATCT'
            """
            if len(seq) <= 6:
                return seq
            else:
                return seq[:6]

        if Counter(headseq(seq))['A'] != 6:
            return False
        return True

    def pAsite(self, strand):
        """
        The read's strand must be same to peak's strand, then return the start site and count
        version alpha: just for single end 10X (R2) peaks detecting
        :param strand:
        :return:
        """
        revornot = self.readline.is_reverse
        if not revornot and strand == '+':

            return self.readline.reference_end + 1

        elif revornot and strand == '-':
            return self.readline.reference_start + 1

        else:
            return False

    def __str__(self):
        return str(self.readline)

    def __getattr__(self, item):
        return getattr(self.readline, item)


def processBarcode(barcodefile):
    """
    For 10X results, the barcode contained the '-1' suffix.
    :param file:
    :return:
    """
    priorinfo = {}
    ccount = 0
    with open(barcodefile) as IN:
        for i in IN.readlines():
            bc = i[:-3]
            ccount += 1
            priorinfo[bc] = ccount
    return priorinfo


def Window(seq, n=2):
    """Return a sliding window over a list
    >>> for i in Window('ATAC', 2):
            print(i)
    'AT'
    'TA'
    'AC'
    """
    it = iter(seq)
    result = list(islice(it, n))

    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + list((elem,))
        yield result


def check_barcode(barcode, pool):
    """Check the barcode whether exsited in known dict, if True then return self
    """
    if barcode in pool:
        return barcode
    else:
        return False


def CheckGap(peaklist):
    """Recieved a peaklist, and check whether intersect bewteen 24 nt
    >>> CheckGap([site, counts])
    """
    sitelength = []
    for i in Window(peaklist, 2):
        sitelength.append(i[1][0] - i[0][0])
    return sitelength


def processBed(bedfile):
    """
    Processing the bed file, and return a list which contains the region information split by chromosomes.
    :param file:
    :return:

    >>> BedParse(bedfile)
    [['chr1\t1\t3'],
    ['chr2\t1\t3'],
    ['chr3\t1\t3']]
    """
    c_chrom = ''
    chrom_list = []
    allinfo = []
    with open(bedfile) as IN:
        for line in IN.readlines():
            line = line.strip()
            chrom = line.split('\t')[0]
            if chrom == c_chrom:
                chrom_list.append(line)
            else:
                if not c_chrom:
                    chrom_list.append(line)
                else:
                    allinfo.extend([chrom_list])
                    chrom_list = []
                    chrom_list.append(line)
                c_chrom = chrom
        allinfo.extend([chrom_list])
    return allinfo


def ContinuousFind(peaklist):
    """Receive a num list and return continuous list in lists.
    >>>ContinuousFind([ 1, 4, 5, 6, 10, 15, 16, 17, 18, 22, 25, 26, 27, 28])
    [[1], [4, 5, 6], [10], [15, 16, 17, 18], [22], [25, 26, 27, 28]]
    """
    finallist = []
    for k, g in groupby(enumerate(peaklist), lambda ix: ix[0] - ix[1]):
        finallist.extend([list(map(itemgetter(1), g))])
    return finallist


def GroupPeak(peaklist, samdic):
    """input: all collapse continuous numbers
        [(site, counts), (site, counts), (site, counts)]
        the sites have been sorted
        samdic was the second collapse sam reads infor
    """

    finalpeak = []

    lastsam = {}

    peakmax = peaklist[-1][0]
    peakmin = peaklist[0][0]
    if peakmax - peakmin >= 24:
        if len(peaklist) == 2:
            return peaklist, samdic

        interval = [peakmin]
        for i in range(ceil((peakmax - peakmin) / 24)):
            interval.append(interval[-1] + 24)

        for subinterval in Window(interval, 2):

            tmp = []
            samtmp = []
            for site in peaklist:
                if subinterval[0] <= site[0] < subinterval[1]:
                    tmp.append(site)
                    samtmp.extend(samdic[site[0]])
                else:
                    if site[0] > subinterval[1]:
                        break
                    else:
                        continue
            if tmp:
                maxsite = max(tmp, key=lambda x: x[1])
                finalpeak.append((maxsite[0], sum([i[1] for i in tmp])))
                lastsam[maxsite[0]] = samtmp
            else:
                continue
        return finalpeak, lastsam
    else:
        if len(peaklist) == 1:
            return peaklist, samdic
        maxsite = max(peaklist, key=lambda x: x[1])

        lastsam[maxsite[0]] = list(chain(*samdic.values()))
        return [(maxsite[0], sum([i[1] for i in peaklist]))], lastsam


def RemoveKey(d, key):
    r = dict(d)
    del r[key]
    return r


def MaxKeys(k, countdic, samdic):
    """Return the max reads count in given regions and the sam reads
    """
    c_key = 0
    c_values = 0
    valuesum = 0
    readsinfo = []
    for i in k:
        if c_key == 0:
            c_key = i
            c_values = countdic[c_key]
            valuesum += countdic[i]
            readsinfo.extend(samdic[i])
        else:
            if c_values > countdic[i]:
                valuesum += countdic[i]
                readsinfo.extend(samdic[i])
                continue
            else:
                c_key = i
                c_values = countdic[i]
                valuesum += countdic[i]
                readsinfo.extend(samdic[i])
    return (c_key, valuesum), {c_key: readsinfo}


def collapsePEAK(peaklist, samdic):
    """input
    """
    dealpeak, dealsamdic = GroupPeak(peaklist, samdic)
    if dealpeak == peaklist:
        sitelength = CheckGap(dealpeak)
        if any(l < 24 for l in sitelength):
            collapse_index = list(map(lambda x: x < 24, sitelength)).index(True)

            tmppeakreads = dealsamdic[dealpeak[collapse_index][0]]
            tmppeakreads += dealsamdic[dealpeak[collapse_index + 1][0]]

            lastsam = RemoveKey(dealsamdic, dealpeak[collapse_index][0])
            lastsam = RemoveKey(lastsam, dealpeak[collapse_index + 1][0])

            maxsite = max(dealpeak[collapse_index:collapse_index + 2], key=lambda x: x[1])
            dealpeak = dealpeak[:collapse_index] + [
                (maxsite[0], sum([i[1] for i in peaklist[collapse_index:collapse_index + 2]]))] + dealpeak[
                                                                                                  collapse_index + 2:]

            lastsam[maxsite[0]] = tmppeakreads
            return collapsePEAK(dealpeak, lastsam)
        return dealpeak, dealsamdic
    else:
        return collapsePEAK(dealpeak, dealsamdic)


def hamming_circle(s, n, alphabet):
    """Generate strings over alphabet whose Hamming distance from s is
    exactly n.

    >>> sorted(hamming_circle('abc', 0, 'abc'))
    ['abc']
    >>> sorted(hamming_circle('abc', 1, 'abc'))
    ['aac', 'aba', 'abb', 'acc', 'bbc', 'cbc']
    >>> sorted(hamming_circle('aaa', 2, 'ab'))
    ['abb', 'bab', 'bba']

    """
    for positions in combinations(range(len(s)), n):
        for replacements in product(range(len(alphabet) - 1), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]
                else:
                    cousin[p] = alphabet[r]
            yield ''.join(cousin)




def processBam(args):
    bamfile, peakinfo, fileindex,DIR= args
    baminfo = pysam.AlignmentFile(bamfile, 'rb', check_sq=False)
    # bamout = pysam.AlignmentFile(os.path.join('pAsupport', '{}.bam'.format((fileindex))), 'wb', header=baminfo.header)
    PAsitecollet = defaultdict(dict)
    chrom, _, _, _, _, _ = peakinfo[0].strip().split('\t')
    pacount = open('{}/{}.txt'.format(DIR,chrom),'w')
    logging.info('Processing the chromosome {}!'.format(chrom))
    print('aaaaa')
    for region in peakinfo:
        # ['chr1\tstart\tend...']
        print(region)
        chrom, st, en, noneinfo, enscore, strand = region.strip().split('\t')
        csite = []
        firstbam = defaultdict(list)
        for read in baminfo.fetch(str(chrom), int(st), int(en), multiple_iterators=True):
            read = ReadParser(read)
            if read.softclip:
                pAsite = read.pAsite(strand)
                try:
                    if int(st) <= pAsite <= int(en):
                        csite.append(pAsite)
                        firstbam[pAsite].append(read)
                    else:
                        continue
                except:
                    pass
        csite = sorted(csite)
        unique_site = sorted(list(set(csite)))
        csite_count = Counter(csite)
        merge_contin = []
        secosam = {}

        for i in ContinuousFind(unique_site):
            if len(i) == 1:
                merge_contin.append((i[0], csite_count[i[0]]))
                secosam[i[0]] = firstbam[i[0]]
            else:
                max_conti, max_conti_reads = MaxKeys(i, csite_count, firstbam)
                merge_contin.append(max_conti)
                secosam.update(max_conti_reads)
        if not merge_contin:
            continue

        finalmerge, finalsam = collapsePEAK(merge_contin, secosam)

        for site, reads in finalsam.items():

            # BCtmp = defaultdict(list)
            cPA = chrom + ':' + str(site) + ':' + strand
            pacount.write('{}\t{}\n'.format(cPA,len(reads)))
            # R1list = []
            # PaLengthSet = []
            # PaAnnotation = []
            # PaLocation = []
            #
            # PaAnnotationPool = []
            # PaLocationPool = []  # if not found in a gene, then return a location. eg. intronic, intergenic

            # for read in reads:
            #     # readBC = read.get_tag('CR')
            #     # cBC = check_barcode(readBC, priorinfo)
            #     # if not cBC: continue
            #     bamout.write(read.readline)

            #     # 把这里的+1直接改成添加barcode的字典，随后计数，这样就完成了对umi的去重复
            #     # 此处的barcode我们使用数字来代替， 10-13
            #     BCtmp[priorinfo[cBC]] = [read.get_tag('UR')]
            # for k, v in BCtmp.items():
            #     v = set(v)
            #     pool = set()  # checking whether ignore these one distance value
            #     collpasedic = {}  # fianl result
            #     for subumi in list(v):
            #         statu = False
            #         if subumi in pool: continue
            #         for onedistance in list(hamming_circle(subumi, 1, "ATGC")):
            #             if onedistance in v:
            #                 statu = True
            #                 if subumi not in collpasedic:
            #                     collpasedic[subumi] = [onedistance]
            #                 else:
            #                     collpasedic[subumi].append(onedistance)
            #                 pool.add(onedistance)
            #             else:
            #                 continue
            #         if not statu:
            #             collpasedic[subumi] = []
            #         else:
            #             continue
            #     umi = len(collpasedic)
            #     PAsitecollet[cPA][k] = umi
    bamout.close()
    pacount.close()
    logging.info('done!')
    return PAsitecollet

def checkdir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def wraperProcess(arg):
    return processBam(arg)


def main(*file):
    bamfile, peakfile, DIR = file[0]
    # barcode = processBarcode(barcodefile)
    peakfile = processBed(peakfile)

    pool = Pool(processes=20)
    res = []
    checkdir(DIR)
    print(len(peakfile))
    for x in range(len(peakfile)):
        arg = [bamfile, peakfile[x], x,DIR]
        res.append(pool.apply_async(wraperProcess, (arg,)))
    pool.close()
    pool.join()

    # WriteInfo(DIR, res)


if __name__ == '__main__':
    main(sys.argv[1:])
