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

from flaimapper.FlaiMapper import FlaiMapper
from flaimapper.CLI import CLI
from flaimapper.Data import *

import unittest


class TestFlaiMapper(unittest.TestCase):
    def test_01_a(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01])
        
        args.parameters.left_padding = 1
        args.parameters.right_padding = 1
        
        flaimapper = FlaiMapper(args.alignment_file,args.verbosity)
        
        i = 0
        for region in flaimapper.regions(args.parameters):
            i += 1
        
        self.assertEqual(i, 3)
    
    def test_01_b(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01])
        
        args.parameters.left_padding = 2
        args.parameters.right_padding = 2
        
        flaimapper = FlaiMapper(args.alignment_file,args.verbosity)
        
        i = 0
        for region in flaimapper.regions(args.parameters):
            i += 1
        
        self.assertEqual(i, 3)
    
    def test_01_c(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01])
        
        args.parameters.left_padding = 3
        args.parameters.right_padding = 3
        
        flaimapper = FlaiMapper(args.alignment_file,args.verbosity)
        
        i = 0
        for region in flaimapper.regions(args.parameters):
            i += 1
        
        self.assertEqual(i, 3)
    
    def test_01_d(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01])
        
        args.parameters.left_padding = 3
        args.parameters.right_padding = 4
        
        flaimapper = FlaiMapper(args.alignment_file,args.verbosity)
        
        i = 0
        for region in flaimapper.regions(args.parameters):
            i += 1
        
        self.assertEqual(i, 2)
    
    def test_01_e(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01])
        
        args.parameters.left_padding = 4
        args.parameters.right_padding = 3
        
        flaimapper = FlaiMapper(args.alignment_file,args.verbosity)
        
        i = 0
        for region in flaimapper.regions(args.parameters):
            i += 1
        
        self.assertEqual(i, 2)
    
    def test_01_f(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01])
        args.parameters.left_padding = 4
        args.parameters.right_padding = 4
        
        flaimapper = FlaiMapper(args.alignment_file,args.verbosity)
        
        i = 0
        for region in flaimapper.regions(args.parameters):
            i += 1
        
        self.assertEqual(i, 2)
    
    
    #def test_10(self):
        #args = CLI([TESTS_EXAMPLE_ALIGNMENT_01])
        
        #flaimapper = FlaiMapper(args.alignment_file,args.verbosity)
        
        ## Run analysis
        #flaimapper.run(args.fasta_handle,args.parameters.matrix)

        #flaimapper.write(args.format,args.output)
        
        #self.assertEqual("status","complete")


def main():
    unittest.main()

if __name__ == '__main__':
    main()
