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



import os,re,random,operator,argparse,sys


from flaimapper.Read import Read
from flaimapper.ncRNAfragment import ncRNAfragment
from flaimapper.MaskedRegion import MaskedRegion



class SSLMParser(MaskedRegion):
	"""parseNcRNA is a class that parses the SSLM alignment files.
	"""
	def parse_reads(self):
		"""parse
		
		----
		"""
		
		previous_line = ""
		
		for filename in self.get_alignment_files():
			i = 0
			with open(filename,'r') as fh:
				for line in fh:
					line = line.strip()
					
					if(i % 2 == 1):
						if(i == 1):
							self.sequence = line
						else:
							# previous_line = ">fasta name _hits etc"
							#          line = "-----ACTG-----"
							
							k = previous_line.lower().find('_hits') 
							
							if(k > -1):
								name = previous_line[1:k]
								numberofhits = int(previous_line[k+5::])
							else:
								name = previous_line[::-1].lstrip(">")
								numberofhits = 1
							
							start_pos = self.get_start_position(line)
							stop_pos = self.get_stop_position(line)
							
							for j in range(numberofhits):
								yield Read(start_pos,stop_pos,name,line[start_pos:stop_pos])
					else:
						previous_line = line
					
					i += 1
	
	def get_alignment_files(self):
		for alignment_directory in self.alignments:
			with open(alignment_directory+"/idreadable.txt",'rU') as fh:
				for line in fh:
					line = line.strip()
					if(line != ""):
						line = line.split("\t")
						if(line[0].lstrip(">") == self.name):
							yield alignment_directory+"/validated/"+line[1]+".fa"

	
	def get_start_position(self,read,extention='-'):
		"""
		Finds start-position of a read according to lines of the following format:
		-----TACCCTGTAGAGCCGAATTTGT-----
		     *
		
		This example will return: 5 because at the 5th position is a 'T' (notice that we count from 0)
		
		----
		@return:
		@rtype: integer
		"""
		
		return len(read)-len(read.lstrip(extention))
	
	def get_stop_position(self,read,extention='-'):
		"""
		Finds stop-position of a read according to lines of the following format:
		-----TACCCTGTAGAGCCGAATTTGT-----
		                           *
		
		This example will return: 27 because at the 27th position is a '-' (notice that we count from 0)
		
		----
		@return:
		@rtype: integer
		"""
		
		return len(read.rstrip(extention))
