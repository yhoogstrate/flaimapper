#!/usr/bin/env python

from flaimapper.miRBase import miRBase
from flaimapper.ncRNA import ncRNA

from flaimapper.utils import fasta_entry_names
from flaimapper.utils import parse_gff
from flaimapper.utils import link_mirbase_to_ncrnadb09
from flaimapper.FlaiMapperObject import FlaiMapperObject

import sys
tmp_dir = sys.argv[1].rstrip("/")+"/"

verbosity = "quiet"

miRNAs = miRBase("../../../../share/annotations/miRBase_20/miRNA.dat")
ncrna_library_names = fasta_entry_names("../../../../share/annotations/ncRNA_annotation/ncrnadb09.fa")
regions = parse_gff("../../../../share/annotations/ncRNA_annotation/ncrnadb09.gtf")
links = link_mirbase_to_ncrnadb09(miRNAs,ncrna_library_names)			# Crosslink miRBase with reference ncRNAs

dataset_id = "SRP041082"

experiments = ['SRR1232072','SRR1232073']
for experiment in experiments:
	alignments = ["../../../../share/small_RNA-seq_alignments/"+dataset_id+"/"+experiment]
	
	# Load flaimapper
	flaimapper = FlaiMapperObject('sslm',verbosity)
	for alignment in alignments:
		flaimapper.add_alignment(alignment)
	results = flaimapper.count_reads_per_region(miRNAs,links,regions,10)
	
	keys = sorted(results.keys())
	keys.append("stacked")
	
	# Export found-backs
	fh = open(tmp_dir+"validation_FlaiMapper_miRBase_"+dataset_id+"_"+experiment+"_sensitivity.txt","w")
	
	factors = ["predicted","not_predicted_with_reads","not_predicted_no_reads"]
	fh.write("factor\t"+"\t".join(keys)+"\n")
	for factor in factors:
		line = factor
		line += "\t"+str(results[keys[0]][factor])
		line += "\t"+str(results[keys[1]][factor])
		line += "\t"+str(results[keys[0]][factor]+results[keys[1]][factor])
		fh.write(line+"\n")
	
	fh.close()
	
	
	# Export Errors
	fh = open(tmp_dir+"validation_FlaiMapper_miRBase_"+dataset_id+"_"+experiment+"_offset.txt","w")
	
	header = "error"
	for k in keys:
		for prime in ["_5p","_3p"]:
			header += "\t"+k+prime
	
	fh.write(header+"\n")
	xkeys = ["<-5",-5,-4,-3,-2,-1,0,1,2,3,4,5,">5"]
	
	for n in xkeys:
		line = str(n)
		
		for k in keys:
			for prime in ["_5p","_3p"]:
				if(k != "stacked"):
					line += "\t"+str(results[k]["error"+prime][n])
				else:
					line += "\t"+str(results[keys[0]]["error"+prime][n]+results[keys[1]]["error"+prime][n])
		
		fh.write(line+"\n")
	fh.close()
