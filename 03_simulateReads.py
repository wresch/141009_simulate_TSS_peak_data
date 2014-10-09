#! /usr/bin/env python

import sys
from Bio import SeqIO
import numpy

fq = "@{sid:05d}\n{seq}\n+\nCCCFFFFFHHHHHJJIJJJJJJJJJIJJJJJJIIJJJJJJJJJJJJJJJJ"

def readGtf(fn):
    tss = []
    with open(fn) as gtf:
        for line in gtf:
            chrom, _, _, starts, ends, _, strand, _ = line.split("\t", 7)
            if strand == "+":
                tss.append((int(starts) - 1, "+"))
            else:
                tss.append((int(ends) - 1, "-"))
    return tss

def readFa(fn):
    return SeqIO.read(fn, "fasta").seq

def sim_reads(tss_list, frag_size, fe, peakwidth, libsize, bpseq):
    hfs = frag_size // 2
    ntss = len(tss_list)
    slen = len(bpseq)
    
    nptss = libsize // ntss
    if nptss % 2 == 1:
        nptss += 1
    hptss = nptss // 2
    spec_reads  = nptss // (fe + 1) * fe  #number of specific reads per tss
    uspec_reads = nptss - spec_reads # number of non-specific reads per tss

    shift = (numpy.random.poisson(50, nptss) - 50 + frag_size) // 2  # half frag size for shifting reads
    peak_var = (peakwidth // 2) ** 2

    n = 0
    for pos, strand in tss_list:
        frag_centers = numpy.concatenate([numpy.random.randint(low = -2000, high = 2000, size = uspec_reads) + pos,
                                         numpy.random.poisson(peak_var, spec_reads) - peak_var + pos])
        # plus strand reads
        plusi = numpy.random.randint(low = 0, high = nptss - 1, size = hptss)
        for start_5p in frag_centers[plusi] - shift[plusi]:
            rleft = start_5p
            rright = start_5p + 50
            if rright < slen and rleft >= 0:
                n += 1
                print fq.format(sid = n, seq = bpseq[rleft:rright])
        minusi = numpy.random.randint(low = 0, high = nptss - 1, size = hptss)
        for start_5p in frag_centers[minusi] + shift[minusi]:
            rleft = start_5p - 49
            rright = start_5p + 1
            if rright < slen and rleft >= 0:
                n += 1
                print fq.format(sid = n, seq = bpseq[(start_5p - 49):(start_5p + 1)].reverse_complement())
         


gtf_file = sys.argv[1]
fa_file  = sys.argv[2]
gtf = readGtf(gtf_file)
fa  = readFa(fa_file)

sim_reads(tss_list  = gtf,
          frag_size = 150,
          fe        = 10,
          peakwidth = 100,
          libsize   = 1000000,
          bpseq     = fa)
