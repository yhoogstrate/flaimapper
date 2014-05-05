#!/usr/bin/env python

import sys
import fileinput
import re

sys.path.append("../../../../src")

from miRBase import *

from FragmentContainer import FragmentContainer
from FragmentFinder import FragmentFinder
from AlignmentParser import AlignmentParser
from AlignmentDirectory import AlignmentDirectory
from FlaiMapperObject import FlaiMapperObject



miRNAs = miRBase("../../../../share/annotations/miRBase_20/miRNA.dat")
verbosity = "quiet"

for sd in [0.05]:#[0.05,0.1,0.2,0.35,0.5,1,1.5,2] <- if every settings should be plotted
	for dist in [26]:#range(1,51) <- if every settings should be plotted
		param = 'min-dist_'+str(dist)+"_sd_"+str(sd)
		print "processing "+param
		filename = "../../../../output/BlockBuster/SRP002175/"+param+".txt"
		
		blockbuster_clusters = {}
		blockbuster_blocks = {}
		
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
		results = flaimapper.count_reads_per_region_custom_table(miRNAs,links,flaimapper.sequences,10)
		
		
		keys = sorted(results.keys())
		keys.append("stacked")
		
		# Export found-backs
		fh = open("1_bottom_validation_BlockBuster_"+param+"_miRBase_SRP00217501_sensitivity.txt","w")
		
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
		fh = open("1_bottom_validation_BlockBuster_"+param+"_miRBase_SRP00217501_offset.txt","w")
		
		header = "error"
		for k in keys:
			for prime in ["_5p","_3p"]:
				header += "\t"+k+prime
		
		fh.write(header+"\n")
		for n in ["<-5",-5,-4,-3,-2,-1,0,1,2,3,4,5,">5"]:
			line = str(n)
			
			for k in keys:
				for prime in ["_5p","_3p"]:
					if(k != "stacked"):
						line += "\t"+str(results[k]["error"+prime][n])
					else:
						line += "\t"+str(results[keys[0]]["error"+prime][n]+results[keys[1]]["error"+prime][n])
			
			fh.write(line+"\n")
		
		fh.close()
		
		del(flaimapper,blockbuster_clusters,predicted_fragments)
