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


from flaimapper.CLI import CLI_sslm2sam
from flaimapper.FragmentContainer import FragmentContainer

import unittest,subprocess,os


class TestFlaiMapper(unittest.TestCase):
    def test_01_a(self):
        os.chdir("../share/small_RNA-seq_alignments/SRP006788/")
        subprocess.call(["tar","-xzf","SRR207111_HeLa18-30.tar.gz"])
        
        args = CLI_sslm2sam(['-o','test.sam','SRR207111_HeLa18-30'])
        sslm2bed_converter = FragmentContainer()
        sslm2bed_converter.convert_to_sam(args.output)
"""
    sslm2bed_converter = FlaiMapperObject('sslm',args.verbosity)
    for alignment_directory in args.alignment_directories:
        sslm2bed_converter.add_alignment(alignment_directory)
    
    regions = parse_gff(args.mask)
    
    sslm2bed_converter.convert_to_sam(regions,args.output)

        """



def main():
    unittest.main()

if __name__ == '__main__':
    main()
