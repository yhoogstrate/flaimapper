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

import os,tempfile

def fasta_entry_names(fasta_file):
	names = {}
	with open(fasta_file) as fh:
		for line in fh:
			line = line.strip()
			if(line != '' and line[0] == '>'):
				names[line[1:]] = True
	return names.keys()

def parse_gff(gff_file):
	"""2015-mar-20: Removed the Tabix library because of incompatibility
	issues.
	"""
	
	regions = []
	
	with open(gff_file,'r') as fh:
		for line in fh:
			line = line.strip()
			if(len(line) > 0 and line[0] != '#'):
				region = line.split('\t')
				regions.append((region[0],int(region[3]),int(region[4])-1,0))
	
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
