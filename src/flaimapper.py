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
from FragmentFinder import FragmentFinder
from AlignmentParser import AlignmentParser
from AlignmentDirectory import AlignmentDirectory
from FlaiMapperObject import FlaiMapperObject



if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	
	group = parser.add_mutually_exclusive_group()
	group.add_argument("-v", "--verbose", action="store_true")
	group.add_argument("-q", "--quiet", action="store_true")
	
	parser.add_argument("-o","--output",help="output filename; '-' for stdout",default="-")
	parser.add_argument("-f","--format",help="file format of the output: [1: table; per fragment], [2: table; per ncRNA], [3: genbank]",type=int,default=1)
	
	parser.add_argument("alignment_directories",nargs='+')
	
	args = parser.parse_args()
	args.verbosity = "normal"
	if(args.verbose == "verbose"):
		args.verbosity = "verbose"
	elif(args.quiet):
		args.verbosity = "quiet"
	
	
	
	flaimapper = FlaiMapperObject(args.verbosity)
	for alignment_directory in args.alignment_directories:
		flaimapper.add_alignment_directory(AlignmentDirectory(alignment_directory,args.verbosity))
	
	flaimapper.run()
	flaimapper.write(args.format,args.output)
