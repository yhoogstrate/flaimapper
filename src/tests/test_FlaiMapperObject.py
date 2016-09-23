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

from flaimapper.FlaiMapperObject import FlaiMapperObject
from flaimapper.CLI import CLI

import unittest,pysam


class TestCLI(unittest.TestCase):
    def test_01(self):
        args = CLI(["test-file.bam"])
        
        # Load BAM Files
        flaimapper = FlaiMapperObject('bam',args.verbosity)
        for alignment_file in args.alignment_files:
            flaimapper.add_alignment(alignment_file)
        
        # The genomic regions of the precursor sequence(s).
        fasta_ref = pysam.Fastafile(args.fasta)
        
        # Run analysis
        flaimapper.run(regions,fasta_ref,args.parameters.matrix)
        flaimapper.write(args.format,args.output)
        
        
        self.assertEqual("status","complete")


def main():
    unittest.main()

if __name__ == '__main__':
    main()
