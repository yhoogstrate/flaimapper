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



from SSLM2BED import SSLM2BED
from AlignmentDirectory import AlignmentDirectory



if __name__ == "__main__":
	"""
	This program converts the alignments of the used format in the
	article to the BED format.
	"""
	parser = argparse.ArgumentParser()
	
	group = parser.add_mutually_exclusive_group()
	group.add_argument("-v", "--verbose", action="store_true")
	group.add_argument("-q", "--quiet", action="store_true")
	parser.add_argument("-o","--output",help="output BED-filename; '-' for stdout",default="-")
	parser.add_argument("alignment_directories",nargs='+')
	
	args = parser.parse_args()
	args.verbosity = "normal"
	if(args.verbose):
		args.verbosity = "verbose"
	elif(args.quiet):
		args.verbosity = "quiet"
	
	
	sslm2bed_converter = SSLM2BED(args.verbosity)
	for alignment_directory in args.alignment_directories:
		sslm2bed_converter.add_alignment_directory(AlignmentDirectory(alignment_directory,args.verbosity))
	
	sslm2bed_converter.convert(args.output)
