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

import flaimapper
from flaimapper.CLI import CLI
from flaimapper.FlaiMapper import FlaiMapper
from flaimapper.Data import *


import unittest,logging
logging.basicConfig(format=flaimapper.__log_format__, level=logging.DEBUG)



class TestFunctional(unittest.TestCase):
    def test_01_a(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,'--verbose'])
        
        fm = FlaiMapper(args.alignment_file)
        
        # First check whether it's looking in the right regions
        i = 0
        for region in fm.regions(args.parameters):
            if i == 0:
                self.assertEqual(region[0], 'chr1')
                self.assertEqual(region[1], max(0, 13-1-args.parameters.left_padding))
                self.assertEqual(region[2], 43+20+args.parameters.right_padding)
            elif i == 1:
                self.assertEqual(region[0], 'chr2')
                self.assertEqual(region[1], 54-1-args.parameters.left_padding)
                self.assertEqual(region[2], 54+26+args.parameters.right_padding)
            else:
                self.asserTrue(False,"Race condition?")
            i += 1
        self.assertEqual(i, 2)
        
        # Then detect the fragments
        i = 0
        fm.run(args.fasta_handle,args.parameters)
        for fragment in fm:
            if i == 0:
                self.assertEqual(fragment.masked_region[0], 'chr1')
                self.assertEqual(fragment.start, 13)
                self.assertEqual(fragment.stop, 36)
            elif i == 1:
                self.assertEqual(fragment.masked_region[0], 'chr1')
                self.assertEqual(fragment.start, 43)
                self.assertEqual(fragment.stop, 62)
            elif i == 2:
                self.assertEqual(fragment.masked_region[0], 'chr2')
                self.assertEqual(fragment.start, 54)
                self.assertEqual(fragment.stop, 79)
            else:
                self.asserTrue(False,"Race condition?")
            i += 1
        self.assertEqual(i, 3)
    
    def test_01_b(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,'-p',TESTS_FUNCTIONAL_DUCK7_PARAMS,'--verbose'])
        
        fm = FlaiMapper(args.alignment_file)
        
        # First check whether it's looking in the right regions
        i = 0
        for region in fm.regions(args.parameters):
            if i == 0:
                self.assertEqual(region[0], 'chr1')
                
                self.assertEqual(region[1], 13-1-args.parameters.left_padding)
                self.assertEqual(region[2], 43+20+args.parameters.right_padding)
            elif i == 1:
                self.assertEqual(region[0], 'chr2')
                
                self.assertEqual(region[1], 54-1-args.parameters.left_padding)
                self.assertEqual(region[2], 54+26+args.parameters.right_padding)
            else:
                self.asserTrue(False,"Race condition?")
            i += 1
        
        self.assertEqual(i, 2)
        
        # Then detect the fragments
        i = 0
        fm.run(args.fasta_handle,args.parameters)
        for fragment in fm:
            if i == 0:
                self.assertEqual(fragment.masked_region[0], 'chr1')
                self.assertEqual(fragment.start, 13)
                self.assertEqual(fragment.stop, 36)
            elif i == 1:
                self.assertEqual(fragment.masked_region[0], 'chr1')
                self.assertEqual(fragment.start, 43)
                self.assertEqual(fragment.stop, 62)
            elif i == 2:
                self.assertEqual(fragment.masked_region[0], 'chr2')
                self.assertEqual(fragment.start, 54)
                self.assertEqual(fragment.stop, 79)
            else:
                self.asserTrue(False,"Race condition?")
            i += 1
        self.assertEqual(i, 3)

    def test_01_c(self):
        """
        This should filter out the second fragment because it's 26bp away.
        """
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,'-p',TESTS_FUNCTIONAL_DUCK26_PARAMS,'--verbose'])
        
        fm = FlaiMapper(args.alignment_file)
        
        # First check whether it's looking in the right regions
        i = 0
        for region in fm.regions(args.parameters):
            if i == 0:
                self.assertEqual(region[0], 'chr1')
                
                self.assertEqual(region[1], max(0, 13-1-args.parameters.left_padding))
                self.assertEqual(region[2], 43+20+args.parameters.right_padding)
            elif i == 1:
                self.assertEqual(region[0], 'chr2')
                
                self.assertEqual(region[1], 54-1-args.parameters.left_padding)
                self.assertEqual(region[2], 54+26+args.parameters.right_padding)
            else:
                self.asserTrue(False,"Race condition?")
            i += 1
        self.assertEqual(i, 2)
        
        # Then detect the fragments
        i = 0
        fm.run(args.fasta_handle,args.parameters)
        for fragment in fm:
            if i == 0:
                self.assertEqual(fragment.masked_region[0], 'chr1')
                self.assertEqual(fragment.start, 13)
                self.assertEqual(fragment.stop, 36)
            elif i == 1:
                self.assertEqual(fragment.masked_region[0], 'chr2')
                self.assertEqual(fragment.start, 54)
                self.assertEqual(fragment.stop, 79)
            else:
                self.asserTrue(False,"Race condition?")
            i += 1
        self.assertEqual(i, 2)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
