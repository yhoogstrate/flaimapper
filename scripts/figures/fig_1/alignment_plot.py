#!/usr/bin/env python

alignment_file = "../../../share/small_RNA-seq_alignments/SRP006788/SRR207111_HeLa18-30/validated/file237.fa"
output_table = "alignment_plot.txt"

fh = open(alignment_file,"r")
lines = fh.read().split("\r")
fh.close()

info = lines[0]
sequence = lines[1]

nlines = len(lines) - 2
nlines = nlines - (nlines % 2)

print "Doing sequence: "+info
print sequence
print "  of length: "+str(len(sequence))

def parse_readcount(name):
	return int(name.lower()[::-1].split("_hits"[::-1])[0][::-1])

def parse_start_pos(alignment,offset=1):# Here counting starts at zero, in R it starts at 1
	return len(alignment) - len(alignment.lstrip('-')) + offset

def parse_stop_pos(alignment,offset=1):# Here counting starts at zero, in R it starts at 1
	return len(alignment.rstrip('-'))+offset

fh = open(output_table,"w")
fh.write("hits\tstart\tstop\n")
for i in range(nlines/2):
	j = (i * 2) + 2
	
	seq_name = lines[j].strip()
	seq_alig = lines[j+1].strip()
	
	fh.write(str(parse_readcount(seq_name))+"\t"+str(parse_start_pos(seq_alig))+"\t"+str(parse_stop_pos(seq_alig))+"\n")

fh.close()
