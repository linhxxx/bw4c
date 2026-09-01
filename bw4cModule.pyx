#!/usr/bin/env python


from __future__ import print_function
import tempfile
import sys
import os
import os.path as op
from subprocess import Popen, PIPE
import argparse
from subprocess import check_call
from operator import itemgetter
from itertools import groupby, repeat, chain, islice
import re
import pysam
from multiprocessing import Process

__version__ = "0.3"

try:
    from itertools import izip
    import string
    maketrans = string.maketrans
except ImportError: # python3
    izip = zip
    maketrans = str.maketrans
import toolshed
from toolshed import nopen, reader, is_newer_b

def comp(s, _comp=maketrans('ATCG', 'TAGC')):
        return s.translate(_comp)

def calASse(bamFile):
    bamIter = pysam.AlignmentFile(bamFile, "rb") 
    asOT = os.popen("zstd -c > "  + bamFile + ".AS.zst", "w")

    for aln in bamIter:
        if not aln.is_supplementary and not aln.is_secondary: 
            asOT.write("{0}\t{1}\n".format(aln.query_name, aln.get_tag("AS"))) 

    bamIter.close()
    asOT.close()

def calASpe(bamFile):
    bamIter = pysam.AlignmentFile(bamFile, "rb") 
    asOT = os.popen("zstd -c > "  + bamFile + ".AS.zst", "w")
    beforeReadname = ""
    beforeScore = 0
    for aln in bamIter:
        if not aln.is_supplementary and not aln.is_secondary: 
            if aln.query_name != beforeReadname: # first meet, just save
                beforeReadname = aln.query_name
                beforeScore = int(aln.get_tag("AS"))
            else:
#                    assert (beforeReadname == aln.query_name, "not same read name in line BAM file")  # no need always true
                    # lenght of queryname is not fixed for different data, not easy to pack into binary, so only use plain txt
                asOT.write("{0}\t{1}\n".format(aln.query_name, int(aln.get_tag("AS")) + beforeScore)) 
    bamIter.close()
    asOT.close()


def alnSelection(bamIter, indexIter, bamOT):
    index_name = (indexIter.readline()).rstrip()
    last_equal_name = ""
    for aln in bamIter:
        if aln.query_name == index_name: # first meet
            bamOT.write(aln)
            last_equal_name = index_name
            index_name = (indexIter.readline()).rstrip() # idexIter always true although all data is readed, but readline return ""
        else:
            if aln.query_name == last_equal_name:
                bamOT.write(aln)


def convert_reads(fq1s, fq2s, out=sys.stdout):

    for fq1, fq2 in zip(fq1s.split(","), fq2s.split(",")):
        sys.stderr.write("converting reads in %s,%s\n" % (fq1, fq2))
        fq1 = nopen(fq1)

        #examines first five lines to detect if this is an interleaved fastq file
        first_five = list(islice(fq1, 5))

        r1_header = first_five[0]
        r2_header = first_five[-1]

        if r1_header.split(' ')[0] == r2_header.split(' ')[0]:
            already_interleaved = True
        else:
            already_interleaved = False

        q1_iter = izip(*[chain.from_iterable([first_five,fq1])] * 4)

        if fq2 != "NA":
            fq2 = nopen(fq2)
            q2_iter = izip(*[fq2] * 4)
        else:
            if already_interleaved:
                sys.stderr.write("detected interleaved fastq\n")
            else:
                sys.stderr.write("WARNING: running bwameth in single-end mode\n")
            q2_iter = repeat((None, None, None, None))

        lt80 = 0

        if already_interleaved:
            selected_iter = q1_iter
        else:
            selected_iter = chain.from_iterable(izip(q1_iter, q2_iter))

        for read_i, (name, seq, _, qual) in enumerate(selected_iter):
            if name is None: continue
            convert_and_write_read(name,seq,qual,read_i%2,out)
            if len(seq) < 80:
                lt80 += 1

    out.flush()
    if lt80 > 50:
        sys.stderr.write("WARNING: %i reads with length < 50\n" % lt80)
        sys.stderr.write("       : this program is designed for long reads\n")
    return 0

def convert_and_write_read(name,seq,qual,read_i,out):

    name = name.rstrip("\r\n").split(" ")[0]
    if name[0] != "@":
        sys.stderr.write("""ERROR!!!!
    ERROR!!! FASTQ conversion failed
    ERROR!!! expecting FASTQ 4-tuples, but found a record %s that doesn't start with "@"
    """ % name)
        sys.exit(1)
    if name.endswith(("_R1", "_R2")):
        name = name[:-3]
    elif name.endswith(("/1", "/2")):
        name = name[:-2]

#    seq = seq.upper().rstrip('\n')
# wont meet lower case seq
    seq = seq.rstrip('\n')
    #char_a, char_b = ['CT', 'GA'][read_i] # 3-letter model
#    char_a, char_b = [['CG', 'TG'], ['CG', 'CA']][read_i] # 4-letter model
# read_i in 0 or 1
    if read_i == 0:
    # keep original sequence as name.
        name = "{0} YS:Z:{1}\tYC:Z:CGTG\n".format(name, seq)
        seq = seq.replace('CG', 'TG')
    else:
        name = "{0} YS:Z:{1}\tYC:Z:CGCA\n".format(name, seq)
        seq = seq.replace('CG', 'CA')

    #sys.stderr.write("test:%s-%s %s\n"%(char_a,char_b,name))
    #sys.stderr.write("changeFQ:%s\n" % seq)
    out.write("".join((name, seq, "\n+\n", qual)))

def as_bam(samfile):
    """
    pfile: either a file or a |process to generate sam output
    fa: the reference fasta
    set_as_failed: None, 'f', or 'r'. If 'f'. Reads mapping to that strand
                      are given the sam flag of a failed QC alignment (0x200).
    """
    #sam_iter = nopen_keep_parent_stdin(pfile, 'r')

    sam_iter = open(samfile, "r")
    for line in sam_iter:
        if not line[0] == "@": break
        handle_header(line)
    else:
        sys.stderr.flush()
        raise Exception("bad or empty fastqs")
    # this line combine the line of the first alignment and the sam_iter
    # BWA output is sort by name already? seems very confident.
    # here groupby + itemgetter(0) is groupby read_name. pair_list  is alignment Bam object
    sam_iter2 = (x.rstrip().split("\t") for x in chain([line], sam_iter))
    for read_name, pair_list in groupby(sam_iter2, itemgetter(0)):
        pair_list = [Bam(toks) for toks in pair_list]
#        print ("[groupby] " +  pair_list)
        for aln in handle_reads(pair_list):
            sys.stdout.write(str(aln) + '\n')

def bam_merge(bamWatsonFile, bamCrickFile, bamOut, pe_or_se):
    # merge watson and crick bam, selection algnment by AS


    #bamW_iter = pysam.AlignmentFile(bamWatsonFile, "rb") # watson
    #bamC_iter = pysam.AlignmentFile(bamCrickFile, "rb") # crick
#    bamW_AS = os.popen("zstd -c > " + bamWatsonFile + ".AS.zst", "w")
#    bamC_AS = os.popen("zstd -c > " + bamCrickFile  + ".AS.zst", "w")

# step1 check AS and judge one 
# the order for primary alignemnt should be the same watson/crick bam, check in bwa output bam already
    if  pe_or_se == "se":
        p1 = Process(target=calASse,args=(bamWatsonFile,))
        p2 = Process(target=calASse,args=(bamCrickFile,))
        p1.start()
        p2.start()
        p1.join()
        p2.join()
    elif pe_or_se == "pe":
        p1 = Process(target=calASpe,args=(bamWatsonFile,))
        p2 = Process(target=calASpe,args=(bamCrickFile,))
        p1.start()
        p2.start()
        p1.join()
        p2.join()

# step Two, slect bigger one Index
    bamW_AS_open = os.popen("zstd -dc " + bamWatsonFile + ".AS.zst", "r")
    bamC_AS_open = os.popen("zstd -dc " + bamCrickFile + ".AS.zst", "r")
    bamW_AS_index = os.popen("zstd -c  > " + bamWatsonFile + ".idx.zst", "w")
    bamC_AS_index = os.popen("zstd -c  > " + bamCrickFile + ".idx.zst", "w")

    for buffW in bamW_AS_open:
        buffC  = bamC_AS_open.readline()
        w_arr = buffW.split("\t")
        c_arr = buffC.split("\t")
        assert (w_arr[0] == c_arr[0], "not same read name in line AS file")
        if int(w_arr[1]) >= int(c_arr[1]):
            bamW_AS_index.write(w_arr[0] + "\n")
        else:
            bamC_AS_index.write(c_arr[0] + "\n")
    bamW_AS_open.close()
    bamC_AS_open.close()
    bamW_AS_index.close()
    bamC_AS_index.close()

# step Three, write bam
    bamW_iter2 = pysam.AlignmentFile(bamWatsonFile, "rb") # watson
    bamOT = pysam.AlignmentFile(bamOut, "wb", template=bamW_iter2) 

    bamW_index_open = os.popen("zstd -dc " + bamWatsonFile + ".idx.zst", "r")
    alnSelection(bamW_iter2, bamW_index_open, bamOT)
    bamW_index_open.close()


    bamC_iter2 = pysam.AlignmentFile(bamCrickFile, "rb") # crick
    bamC_index_open = os.popen("zstd -dc " + bamCrickFile + ".idx.zst", "r")

    alnSelection(bamC_iter2, bamC_index_open, bamOT)
    bamC_index_open.close()

    bamOT.close()
    # step four, clean dir
    os.remove(bamWatsonFile + ".AS.zst")
    os.remove(bamCrickFile + ".AS.zst")
    os.remove(bamWatsonFile + ".idx.zst")
    os.remove(bamCrickFile + ".idx.zst")



def handle_header(line, out=sys.stdout):
    toks = line.rstrip().split("\t")
    if toks[0].startswith("@SQ"):
        sq, sn, ln = toks  # @SQ    SN:fchr11    LN:122082543
        # we have f and r, only print out f
        chrom = sn.split(":")[1]
#        if chrom.startswith('r'): return # if wat and crk seperate align, it will meet ref name begin with "r", and it wont meet duplicate ref name no need to skip it.
        chrom = chrom[1:]
        toks = ["%s\tSN:%s\t%s" % (sq, chrom, ln)]
    if toks[0].startswith("@PG"):
        #out.write("\t".join(toks) + "\n")
        toks = ["@PG\tID:bw4\tPN:bw4\tVN:%s\tCL:\"%s\"" % (
                         __version__,
                         " ".join(x.replace("\t", "\\t") for x in sys.argv))]
    out.write("\t".join(toks) + "\n")

def handle_reads(alns):

    for aln in alns:
        orig_seq = aln.original_seq
        assert len(aln.seq) == len(aln.qual), aln.read
        # don't need this any more.
        aln.other = [x for x in aln.other if not x.startswith('YS:Z')]


        direction = aln.chrom[0]
        aln.chrom = aln.chrom.lstrip('fr') 

        if not aln.is_mapped():
            aln.seq = orig_seq
            continue

        assert direction in 'fr', (direction, aln)
        aln.other.append('YD:Z:' + direction)

#        if set_as_failed == direction:
#            aln.flag |= 0x200

        # here we have a heuristic that if the longest match is not 44% of the
        # sequence length, we mark it as failed QC and un-pair it. At the end
        # of the loop we set all members of this pair to be unmapped

        # as BWA, no fail flag in BAM
  #      if aln.longest_match() < (len(orig_seq) * 0.44):
  #          aln.flag |= 0x200  # fail qc
  #          aln.flag &= (~0x2) # un-pair
  #          aln.mapq = min(int(aln.mapq), 1)

        mate_direction = aln.chrom_mate[0]
        if mate_direction not in "*=":
            aln.chrom_mate = aln.chrom_mate[1:]

        # adjust the original seq to the cigar
        l, r = aln.left_shift(), aln.right_shift()
        if aln.is_plus_read():
            aln.seq = orig_seq[l:r]
        else:
            aln.seq = comp(orig_seq[::-1][l:r])
# same as flag fail mark 
#    if any(aln.flag & 0x200 for aln in alns):
#        for aln in alns:
#            aln.flag |= 0x200
#            aln.flag &= (~0x2)
    return alns

class Bam(object):
    __slots__ = 'read flag chrom pos mapq cigar chrom_mate pos_mate tlen \
            seq qual other'.split()
    def __init__(self, args):
        for a, v in zip(self.__slots__[:11], args):
            setattr(self, a, v)
        self.other = args[11:]
        self.flag = int(self.flag)
        self.pos = int(self.pos)
        self.tlen = int(float(self.tlen))

    def __repr__(self):
        return "Bam({chr}:{start}:{read}".format(chr=self.chrom,
                                                 start=self.pos,
                                                 read=self.read)

    def __str__(self):
        return "\t".join(str(getattr(self, s)) for s in self.__slots__[:11]) \
                         + "\t" + "\t".join(self.other)

    def is_first_read(self):
        return bool(self.flag & 0x40)

    def is_second_read(self):
        return bool(self.flag & 0x80)

    def is_plus_read(self):
        return not (self.flag & 0x10)

    def is_minus_read(self):
        return bool(self.flag & 0x10)

    def is_mapped(self):
        return not (self.flag & 0x4)

    def cigs(self):
        if self.cigar == "*":
            yield (0, None)
            raise StopIteration
        cig_iter = groupby(self.cigar, lambda c: c.isdigit())
        for g, n in cig_iter:
            yield int("".join(n)), "".join(next(cig_iter)[1])

    def cig_len(self):
        return sum(c[0] for c in self.cigs() if c[1] in
                   ("M", "D", "N", "EQ", "X", "P"))

    def left_shift(self):
        left = 0
        for n, cig in self.cigs():
            if cig == "M": break
            if cig == "H":
                left += n
        return left

    def right_shift(self):
        right = 0
        for n, cig in reversed(list(self.cigs())):
            if cig == "M": break
            if cig == "H":
                right += n
        return -right or None

    @property
    def original_seq(self):
        try:
            return next(x for x in self.other if x.startswith("YS:Z:"))[5:]
        except:
            sys.stderr.write(repr(self.other) + "\n")
            sys.stderr.write(self.read + "\n")
            raise

    @property
    def ga_ct(self):
        return [x for x in self.other if x.startswith("YC:Z:")]

    def longest_match(self, patt=re.compile("\d+M")):
        return max(int(x[:-1]) for x in patt.findall(self.cigar))
