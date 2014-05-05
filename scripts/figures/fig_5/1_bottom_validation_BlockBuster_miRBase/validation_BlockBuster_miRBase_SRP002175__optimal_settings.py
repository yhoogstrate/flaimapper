#!/usr/bin/env python

import sys
import fileinput
import re
import os.path

sys.path.append("../../../../src")

from miRBase import *

from FragmentContainer import FragmentContainer
from FragmentFinder import FragmentFinder
from AlignmentParser import AlignmentParser
from AlignmentDirectory import AlignmentDirectory
from FlaiMapperObject import FlaiMapperObject



miRNAs = miRBase("../../../../share/annotations/miRBase_20/miRNA.dat")
verbosity = "quiet"

errors = {}



#for sd in []:
for sd in [0.05,0.1,0.2,0.35,0.5,1,1.5,2]:
	for dist in range(1,51):
		param = 'min-dist_'+str(dist)+"_sd_"+str(sd)
		print "processing "+param
		filename = "../../../../output/BlockBuster/SRP002175/"+param+".txt"
		
		blockbuster_clusters = {}
		blockbuster_blocks = {}
		
		if(not os.path.isfile(filename)):
			print "Could not find BlockBuster output file:\n\t"+filename
		else:
			for line in fileinput.input([filename]):
				line = line.strip()
				if(len(line) > 0):
					if(line[0] == '>'):# Cluster block
						match = re.search(">([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)",line)
						
						blockbuster_cluster = {}
						blockbuster_cluster['start'] = int(match.group(3))
						blockbuster_cluster['stop'] = int(match.group(4))
						blockbuster_cluster['sequence'] = ''
						blockbuster_cluster['expression'] = float(match.group(6))
						
						if(not blockbuster_clusters.has_key(match.group(2))):
							blockbuster_clusters[match.group(2)] = []
						blockbuster_clusters[match.group(2)].append(blockbuster_cluster)
					else:
						match = re.search("([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)",line)
						
						blockbuster_block = {}
						blockbuster_block['start'] = int(match.group(3))
						blockbuster_block['stop'] = int(match.group(4))
						blockbuster_block['sequence'] = ''
						blockbuster_block['expression'] = float(match.group(6))
						
						if(not blockbuster_blocks.has_key(match.group(2))):
							blockbuster_blocks[match.group(2)] = []
						blockbuster_blocks[match.group(2)].append(blockbuster_block)
			
			blockbuster_clusters = blockbuster_blocks
			
			# Crosslink miRBase with reference ncRNAs
			i = 0
			links = {}
			for miRNA in miRNAs.get_miRNAs():
				for name in blockbuster_clusters:
					if(name.lower().find("mir") > -1):
						flaimapper_name_raw = name.split("HUGO-Symbol=")[1].split("&")[0]
						
						if(flaimapper_name_raw in miRNA.get_parameter("aliases")):
							links[name] = miRNA.params["name"]
							
							i += 1
			
			
			# Convert blockbuster into a FlaiMapper object
			flaimapper = FlaiMapperObject(verbosity)
			for ncRNA in blockbuster_clusters:
				predicted_fragments = FragmentFinder(ncRNA,False,False)
				predicted_fragments.results = blockbuster_clusters[ncRNA]
				flaimapper.add_fragments(predicted_fragments)
			
			if(not errors.has_key(sd)):
				errors[sd] = {}
			
			errors[sd][dist] = flaimapper.count_reads_per_region_custom_mse(miRNAs,links,flaimapper.sequences,10)

fh = open("validation_BlockBuster_miRBase__optimal_settings__root_square_error_plateau_start_positions.txt","w")
fh.write("dist")
for sd in errors.keys():
	fh.write("\t"+str(sd))
fh.write("\n")

for dist in errors[sd].keys():
	fh.write(str(dist))
	
	for sd in errors.keys():
		fh.write("\t"+str(errors[sd][dist][0]))
	fh.write("\n")
fh.close()



fh = open("validation_BlockBuster_miRBase__optimal_settings__root_square_error_plateau_stop_positions.txt","w")
fh.write("dist")
for sd in errors.keys():
	fh.write("\t"+str(sd))
fh.write("\n")

for dist in errors[sd].keys():
	fh.write(str(dist))
	
	for sd in errors.keys():
		fh.write("\t"+str(errors[sd][dist][1]))
	fh.write("\n")
fh.close()
