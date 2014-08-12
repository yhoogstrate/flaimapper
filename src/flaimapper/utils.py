#!/usr/bin/env python

"""FlaiMapper: computational annotation of small ncRNA derived fragments using RNA-seq high throughput data

 Here we present Fragment Location Annotation Identification mapper
 (FlaiMapper), a method that extracts and annotates the locations of
 sncRNA-derived RNAs (sncdRNAs). These sncdRNAs are often detected in
 sequencing data and observed as fragments of their  precursor sncRNA.
 Using small RNA-seq read alignments, FlaiMapper is able to annotate
 fragments primarily by peak-detection on the start and  end position
 densities followed by filtering and a reconstruction processes.
 Copyright (C) 2011-2014:
 - Youri Hoogstrate
 - Elena S. Martens-Uzunova
 - Guido Jenster
 
 
 [License: GPL3]
 
 This file is part of flaimapper.
 
 flaimapper is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 flaimapper is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

 Documentation as defined by:
 <http://epydoc.sourceforge.net/manual-fields.html#fields-synonyms>
"""

import os
import pysam

def index_gff(gff_file):
	gff_file = os.path.abspath(gff_file)
	new_idx = gff_file+".tbi"
	if(not os.path.isfile(new_idx)):
		print "making index for gtf/gff file..."
		tmp_name = tempfile.gettempdir()+"/"+os.path.basename(gff_file)
		tmp_name_gz = tmp_name+".gz"
		tmp_name_idx = tmp_name_gz+".tbi"
		
		if(os.path.isfile(tmp_name)):
			os.remove(tmp_name)
		if(os.path.isfile(tmp_name_gz)):
			os.remove(tmp_name_gz)
		if(os.path.isfile(tmp_name_idx)):
			os.remove(tmp_name_idx)
		
		os.symlink(gff_file, tmp_name)
		
		pysam.tabix_index(tmp_name, preset="gff")
		os.rename(tmp_name_idx,new_idx)

def fasta_entry_names(fasta_file):
	names = {}
	with open(fasta_file) as fh:
		for line in fh:
			line = line.strip()
			if(line != '' and line[0] == '>'):
				names[line[1:]] = True
	return names.keys()

def parse_gff(gff_file):
	index_gff(gff_file)
	
	regions = []
	
	for region in pysam.Tabixfile(gff_file).fetch():
		region = region.strip("\r").split("\t")
		regions.append((region[0],int(region[3]),int(region[4])-1,region[6]))
	
	return regions

def link_mirbase_to_ncrnadb09(mirbase,ncrnadb09):
	links = {}
	
	for miRNA in mirbase.get_miRNAs():
		for name in ncrnadb09:
			if(name.lower().find("mir") > -1):
				flaimapper_name_raw = name.split("HUGO-Symbol=")[1].split("&")[0]
				
				if(flaimapper_name_raw in miRNA.get_parameter("aliases")):
					links[name] = miRNA.params["name"]
	
	return links
