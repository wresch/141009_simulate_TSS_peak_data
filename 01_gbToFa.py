#! /usr/bin/env python

import sys
from Bio import SeqIO
from Bio import Entrez

def gbToFa(fn, outfn):
    rec = SeqIO.read(fn, "genbank")
    print "Processing %s [%d nts]" % (rec.id, len(rec))
    rec.id = "chr1"
    SeqIO.write(rec, outfn, "fasta") 


if __name__ == "__main__":
    gb_file = sys.argv[1]
    fa_file = sys.argv[2]
    gbToFa(gb_file, fa_file)
