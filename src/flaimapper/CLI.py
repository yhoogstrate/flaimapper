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


import os,re,operator,argparse,sys,tempfile,textwrap,datetime,logging
import pysam
import flaimapper


def CLI(argv=None):
    from flaimapper.FilterParameters import FilterParameters
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,epilog="Further details can be found in the manual:\n<"+flaimapper.__homepage__+">")
    
    # Writing to stderr, python issue: https://hg.python.org/cpython/rev/ec9a4b77f37b
    parser.add_argument('-V','--version', action='version', version=textwrap.dedent("%(prog)s "+flaimapper.__version__+"\nCopyright (C) 2011-"+str(datetime.datetime.now().year)+" Youri Hoogstrate.\n\nLicense GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\nThis is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n\nFurther details can be found and possible issues can be reported at our homepage:\n<"+flaimapper.__homepage__+">"))
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v","--verbose", action="store_true",default=False)
    group.add_argument("-q","--quiet", action="store_false",default=True)
    
    parser.add_argument("-p","--parameters",required=False,help="File containing the filtering parameters, using default if none is provided")
    
    parser.add_argument("-o","--output",help="output filename; '-' for stdout",default="-")
    parser.add_argument("-f","--format",help="file format of the output: [1: table; per fragment], [2: GTF (default)]",type=int,choices=range(1, 2+1),default=2)
    
    parser.add_argument("-r","--fasta",help="Single reference FASTA file (+faid index) containing all genomic reference sequences")
    
    parser.add_argument("--offset5p",help="Offset in bp added to the exon-type annotations in the GTF file. This offset is used in tools estimating the expression levels (default=4)",type=int,default=4)
    parser.add_argument("--offset3p",help="Offset in bp added to the exon-type annotations in the GTF file. This offset is used in tools estimating the expression levels (default=4)",type=int,default=4)
    
    parser.add_argument("alignment_file",help="indexed SAM or BAM file",nargs=1)
    
    # Parse parameters
    if argv == None:
        args = parser.parse_args()
    else:# Argumented parameters (only for testing)
        args = parser.parse_args(argv)
    
    args.alignment_file = pysam.AlignmentFile(args.alignment_file[0])
    
    if args.fasta != None:
        args.fasta_handle = pysam.Fastafile(args.fasta)
    else:
        args.fasta_handle = None
    
    args.parameters = FilterParameters(args.parameters)
    
    # Set verbosity and logging
    if args.verbose:
        logging.basicConfig(format=flaimapper.__log_format__, level=logging.DEBUG)
        logging.info("Verbose output.")
    elif args.quiet:
        logging.basicConfig(format=flaimapper.__log_format__, level=logging.CRITICAL)
        print "A"
    else:
        logging.basicConfig(format=flaimapper.__log_format__, level=logging.INFO)
        print "B"
    
    return args


def CLI_sslm2sam(argv=None):
    parser = argparse.ArgumentParser()
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v","--verbose", action="store_true",default=False)
    group.add_argument("-q","--quiet", action="store_false",default=True)
    
    parser.add_argument("-o","--output",help="output SAM-filename; '-' for stdout",default="-")
    
    parser.add_argument("sslm_directory",nargs=1,help="SSLM formatted output directories")
    
    # Parse parameters
    if argv == None:
        args = parser.parse_args()
    else:# Argumented parameters (only for testing)
        args = parser.parse_args(argv)
    
    args.sslm_directory = args.sslm_directory[0]
    
    # Set verbosity and logging
    if args.verbose:
        logging.basicConfig(format=flaimapper.__log_format__, level=logging.DEBUG)
        logging.info("Verbose output.")
    elif(args.quiet):
        logging.basicConfig(format=flaimapper.__log_format__, level=logging.CRITICAL)
    else:
        logging.basicConfig(format=flaimapper.__log_format__, level=logging.INFO)
    
    return args
