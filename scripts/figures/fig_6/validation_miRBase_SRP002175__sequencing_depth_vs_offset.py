#!/usr/bin/env python

import sys
sys.path.append("../../../src")



from miRBase import *

from FragmentContainer import FragmentContainer
from FragmentFinder import FragmentFinder
from AlignmentParser import AlignmentParser
from AlignmentDirectory import AlignmentDirectory
from FlaiMapperObject import FlaiMapperObject



verbosity = "quiet"
miRNAs = miRBase("../../../share/annotations/miRBase_20/miRNA.dat")
dataset_id = "SRP002175"

alignments = []
alignments.append(AlignmentDirectory("../../../share/small_RNA-seq_alignments/"+dataset_id+"/SRR038852",verbosity))
alignments.append(AlignmentDirectory("../../../share/small_RNA-seq_alignments/"+dataset_id+"/SRR038853",verbosity))
alignments.append(AlignmentDirectory("../../../share/small_RNA-seq_alignments/"+dataset_id+"/SRR038854",verbosity))
alignments.append(AlignmentDirectory("../../../share/small_RNA-seq_alignments/"+dataset_id+"/SRR038855",verbosity))
alignments.append(AlignmentDirectory("../../../share/small_RNA-seq_alignments/"+dataset_id+"/SRR038856",verbosity))
alignments.append(AlignmentDirectory("../../../share/small_RNA-seq_alignments/"+dataset_id+"/SRR038857",verbosity))
alignments.append(AlignmentDirectory("../../../share/small_RNA-seq_alignments/"+dataset_id+"/SRR038858",verbosity))
alignments.append(AlignmentDirectory("../../../share/small_RNA-seq_alignments/"+dataset_id+"/SRR038859",verbosity))
alignments.append(AlignmentDirectory("../../../share/small_RNA-seq_alignments/"+dataset_id+"/SRR038860",verbosity))
alignments.append(AlignmentDirectory("../../../share/small_RNA-seq_alignments/"+dataset_id+"/SRR038861",verbosity))
alignments.append(AlignmentDirectory("../../../share/small_RNA-seq_alignments/"+dataset_id+"/SRR038862",verbosity))
alignments.append(AlignmentDirectory("../../../share/small_RNA-seq_alignments/"+dataset_id+"/SRR038863",verbosity))

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
results = flaimapper.count_error_with_intensity(miRNAs,links,10)

fh = open("validation_miRBase_SRP002175__sequencing_depth_vs_offset.txt","w")
fh.write("5p_corresponding_reads\t5p_error\t3p_corresponding_reads\t3p_error\n")
for line in results:
	fh.write(str(line['5p'][0])+"\t"+str(line['5p'][1])+"\t"+str(line['3p'][0])+"\t"+str(line['3p'][1])+"\n")
fh.close()

