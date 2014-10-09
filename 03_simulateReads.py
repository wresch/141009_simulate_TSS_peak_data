#! /usr/bin/env python

import sys
import click
from Bio import SeqIO
import numpy

fq = "@{sid:05d}\n{seq}\n+\nCCCFFFFFHHHHHJJIJJJJJJJJJIJJJJJJIIJJJJJJJJJJJJJJJJ"

def readGtf(fh):
    tss = []
    for line in fh:
        chrom, _, _, starts, ends, _, strand, frame, attributes = line.split("\t")
        attrs = (x.strip().split() for x in attributes.split(";") if x.strip() != "")
        attrd = dict((a.strip(), b.strip().strip('"')) for a,b in attrs)
        if attrd["exon_number"] == "1":
            if strand == "+":
                tss.append((chrom, int(starts) - 1, "+"))
            else:
                tss.append((chrom, int(ends) - 1, "-"))
    return tss

def readFa(fh):
    return SeqIO.read(fh, "fasta").seq

def sim_reads(tss_list, frag_size, fe, peakwidth, libsize, bpseqrec):
    bpseq = bpseqrec.seq
    seqid = bpseqrec.id
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
    for chrom, pos, strand in tss_list:
        if chrom != seqid:
            continue
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


@click.command()
@click.argument("gtf_file", type=click.File("rb"))
@click.argument("fa_file", type=click.File("rb"))
@click.option("--frag-size", default=150,
              help="mean ChIP fragment size [150]")
@click.option("--fe", default=10,
              help="""fold enrichment - how many fold more reads originate
from specific sequences than from non-specific sequences [10]""")
@click.option("--peakwidth", default=100,
              help="Width of the chipseq peak [100]")
@click.option("--libsize", default=1000000,
              help="Approximate number of reads to generate [1000000]")
def cli(gtf_file, fa_file, frag_size, fe, peakwidth, libsize):
    """
Given a fasta sequence and a matching GTF file with exons, simulate
ChIP-Seq-like 50 nt fastq data with a fixed quality string from the
sequence. Only one chromosome is used and the chromosome name in the
GTF file must match the first (or only) chromosome in the fasta
file. The exons are used to find TSS sites using the location, strand,
and `exon_number` attribute for each gene."""
    
    tss = readGtf(gtf_file)
    fa  = SeqIO.read(fa_file, "fasta")
    sim_reads(tss_list  = tss,
              frag_size = frag_size,
              fe        = fe,
              peakwidth = peakwidth,
              libsize   = libsize,
              bpseqrec  = fa)


if __name__ == "__main__":
    cli()

