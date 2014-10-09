#! /usr/bin/env python

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# read in genbank record
gb = sys.argv[1]

rec = SeqIO.read(gb, "genbank")

# print GTF format annotation for gosr tssd
prev_tss = 0
n_genes = 0
n_ignored = 0
for feat in rec.features:
    if feat.type != "gene":
        continue
    n_genes += 1
    if feat.strand == 1:
        tss = feat.location.start
        strand = "+"
    elif feat.strand == -1:
        tss = feat.location.end
        strand = "-"
    else:
        sys.exit("bad strand id: %s" % feat.strand)
    if tss - prev_tss < 2000:
        n_ignored += 1
    else:
        gene = feat.qualifiers["gene"][0]
        lt   = feat.qualifiers["locus_tag"][0]
        print "chr1\tgbToGtf\texon\t%d\t%d\t.\t%s\t0\tgene_id \"%s\"; gene_name \"%s\"; exon_number 1; tss_id \"%s\";" % (
            feat.location.start + 1, feat.location.end, strand, lt, gene, gene)
    prev_tss = tss
print >>sys.stderr, "total genes:   %6d" % n_genes
print >>sys.stderr, "ignored genes: %6d" % n_ignored


