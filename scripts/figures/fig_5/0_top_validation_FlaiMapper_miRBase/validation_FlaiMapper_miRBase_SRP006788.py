#!/usr/bin/env python

import sys
sys.path.append("../../../../src")



from miRBase import *

from FragmentContainer import FragmentContainer
from FragmentFinder import FragmentFinder
from AlignmentParser import AlignmentParser
from AlignmentDirectory import AlignmentDirectory
from FlaiMapperObject import FlaiMapperObject



verbosity = "quiet"
miRNAs = miRBase("../../../../share/annotations/miRBase_20/miRNA.dat")
dataset_id = "SRP006788"




experiments = ['SRR207111_HeLa18-30', 'SRR207112_HeLa18-30_RRP40', 'SRR207113_HeLa18-30_AGO1_2', 'SRR207114_HeLa18-30_AGO1_2_RRP40', 'SRR207115_HeLa18-30_XRN1_2', 'SRR207116_HeLa18-30_N']
for experiment in experiments:
	alignments = []
	alignments.append(AlignmentDirectory("../../../../share/small_RNA-seq_alignments/"+dataset_id+"/"+experiment,verbosity))
	
	ncrna_library_names = {}
	for alignment in alignments:
		for name in alignment.index:
			ncrna_library_names[name] = True
	ncrna_library_names = ncrna_library_names.keys()
	
	i = 0
	
	# Crosslink miRBase with reference ncRNAs
	links = {}
	for miRNA in miRNAs.get_miRNAs():
		for name in ncrna_library_names:
			if(name.lower().find("mir") > -1):
				flaimapper_name_raw = name.split("HUGO-Symbol=")[1].split("&")[0]
				
				if(flaimapper_name_raw in miRNA.get_parameter("aliases")):
					links[name] = miRNA.params["name"]
					
					i += 1
	
	# Load flaimapper
	flaimapper = FlaiMapperObject(verbosity)
	for alignment in alignments:
		flaimapper.add_alignment_directory(alignment)
	results = flaimapper.count_reads_per_region(miRNAs,links,10)
	
	
	keys = sorted(results.keys())
	keys.append("stacked")
	
	# Export found-backs
	fh = open("validation_FlaiMapper_miRBase_"+dataset_id+"_"+experiment+"_sensitivity.txt","w")
	
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
	fh = open("validation_FlaiMapper_miRBase_"+dataset_id+"_"+experiment+"_offset.txt","w")
	
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


