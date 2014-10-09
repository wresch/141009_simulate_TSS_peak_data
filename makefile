SHELL := /bin/bash

gb  := nc_000913/nc_000913.gb
fa  := nc_000913/nc_000913.fa
gtf := nc_000913/nc_000913.gtf

$(fa): 01_gbToFa.py $(gb)
	./$^ $@

$(gtf): 02_gbToGtf.py $(gb)
	./$^ > $@

idx := nc_000913/nc_000913
$(idx).bt1.current: $(fa)
	module load bowtie/1.0.0 && \
	  bowtie-build $(fa) $(idx) && \
	  touch $@

coli: $(idx).bt1.current $(gtf)

seq_files:
	mkdir -p seq_files
seq_files/sim.fq.gz: 03_simulateReads.py $(gtf) $(fa)
	./$^ | gzip -c > $@.temp && mv $@.temp $@
sim: seq_files seq_files/sim.fq.gz

aln_files:
	mkdir -p aln_files
aln_files/sim.bam: seq_files/sim.fq.gz $(idx).bt1.current
	module load bowtie/1.0.0 && \
	  gunzip -c $< \
	    | bowtie --phred33-quals --sam --best --strata --all -m1 -n2 --threads=2 -l50 $(idx) - \
	    | samtools view -Sb -F4 - > $@.tmp && mv $@.tmp $@

align: aln_files aln_files/sim.bam


tssd.dat: aln_files/sim.bam $(gtf)
	gosr tssd $^ &> $@.tmp && sed 's/INFO/#INFO/' $@.tmp > $@ && rm $@.tmp