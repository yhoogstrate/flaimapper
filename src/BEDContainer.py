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
 - MSc. Youri Hoogstrate
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



from FragmentContainer import FragmentContainer



class BEDContainer:
	def __init__(self):
		self.genome = {}
	
	def add_chromosome(self,chromosome):
		self.genome[chromosome] = {}
	
	def has_chromosome(self,chromosome):
		return self.genome.has_key(chromosome)
	
	def add_read_to_chromosome(self,position,coverage):
		if(not self.has_chromosome(position['chr'])):
			self.add_chromosome(position['chr'])
		
		for i in range(len(coverage)):
			pos = position['start']+i
			if(not self.genome[position['chr']].has_key(pos)):
				self.genome[position['chr']][pos] = 0
			self.genome[position['chr']][pos] += coverage[i]
	
	def export(self,output = "-"):
		if(output == "-"):
			fh = sys.stdout
		else:
			fh = open(output,"w")
		
		for chromosome in sorted(self.genome.keys()):
			for location in sorted(self.genome[chromosome].keys()):
				## while loop until differencec?
					if(self.genome[chromosome][location] > 0):
						fh.write(chromosome+"\t"+str(location)+"\t"+str(location+1)+"\t"+chromosome+"_"+str(location)+"\t"+str(self.genome[chromosome][location])+"\t-\n")
		fh.close()


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
		
		for ncRNA in self.alignment_directories_indexed.keys():
			if(self.verbosity == "verbose"):
				print "   - Converting: "+ncRNA
			
			read_alignments = AlignmentParser(ncRNA,self.alignment_directories_indexed[ncRNA],self.verbosity)
			self.bed.add_read_to_chromosome({'chr':read_alignments.positions.group(1),'start':int(read_alignments.positions.group(2)),'stop':int(read_alignments.positions.group(3))},read_alignments.coverage)
		
		self.bed.export(output)
