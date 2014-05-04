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


__version_info__ = ('1', '0', '0')
__version__ = '.'.join(__version_info__)
__author__ = 'Youri Hoogstrate'
__homepage__ = 'https://github.com/yhoogstrate/flaimapper'
__license__ = 'GPL3'



import os,re,random,operator,argparse,sys



from BEDContainer import BEDContainer
from FragmentContainer import FragmentContainer
from AlignmentDirectory import AlignmentDirectory
from AlignmentParser import AlignmentParser



class SSLM2BED(FragmentContainer):
	def __init__(self,verbosity):
		self.verbosity = verbosity
		
		self.alignment_directories = []
		self.alignment_directories_indexed = {}
		
		self.bed = BEDContainer()
		
		if(self.verbosity == "verbose"):
			print " - Initiated FlaiMapper Object"
	
	def add_alignment_directory(self,alignment_directory):
		if(alignment_directory.__class__.__name__ != AlignmentDirectory.__name__):
			raise TypeError, "alignment_directory is not of class type " + AlignmentDirectory.__name__ + " but of class type: " + alignment_directory.__class__.__name__
		else:
			self.alignment_directories.append(alignment_directory)
	
	def index(self):
		if(self.verbosity == "verbose"):
			print " - Indexing alignment files"
		
		for alignment_directory in self.alignment_directories:
			if(self.verbosity == "verbose"):
				print "   - Indexing "+alignment_directory.path
			for ncRNA in alignment_directory.index.keys():
				if(not self.alignment_directories_indexed.has_key(ncRNA)):
					self.alignment_directories_indexed[ncRNA] = []
				self.alignment_directories_indexed[ncRNA].append(alignment_directory.index[ncRNA])
	
	def convert(self,output):
		self.index()
		
		if(self.verbosity == "verbose"):
			print " - Running fragment detection"
		
		fh = open(output,"w")
		
		for ncRNA in self.alignment_directories_indexed.keys():
			if(self.verbosity == "verbose"):
				print "   - Converting: "+ncRNA
			
			read_alignments = AlignmentParser(ncRNA,self.alignment_directories_indexed[ncRNA],self.verbosity)
			read_alignments.parse_stats()
			
			reads = read_alignments.parse_reads()
			if(reads):
				for read_id,read in reads.items():
					fh.write(ncRNA+"\t"+str(read['startPos'])+"\t"+str(read['stopPos'])+"\t"+read_id+"\t"+str(read['hits'])+"\t-\n")
		
		fh.close()
