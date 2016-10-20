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
import unittest,filecmp,os,logging,subprocess

logging.basicConfig(format=flaimapper.__log_format__, level=logging.DEBUG)

from flaimapper.FlaiMapper import FlaiMapper
from flaimapper.CLI import CLI
from flaimapper.Data import *


class TestFlaiMapper(unittest.TestCase):
    def test_01_a(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,'--verbose'])
        
        args.parameters.left_padding = 1
        args.parameters.right_padding = 1
        
        flaimapper = FlaiMapper(args)
        
        i = 0
        for region in flaimapper.regions():
            if i == 0:
                self.assertEqual(region.region[0], 'chr1')
                
                self.assertEqual(region.region[1], 13-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 13+20+4+args.parameters.right_padding)
            elif i == 1:
                self.assertEqual(region.region[0], 'chr1')
                
                self.assertEqual(region.region[1], 43-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 43+20+args.parameters.right_padding)
            elif i == 2:
                self.assertEqual(region.region[0], 'chr2')
                
                self.assertEqual(region.region[1], 54-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 54+26+args.parameters.right_padding)
            else:
                self.asserTrue(False,"Race condition?")
            
            i += 1
        
        self.assertEqual(i, 3)

    def test_01_b(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,'-p',TESTS_FUNCTIONAL_DUCK7_PARAMS,'--verbose'])
        
        # Use custom file that must be set to a padding of 7
        self.assertEqual(args.parameters.left_padding, 7)
        self.assertEqual(args.parameters.right_padding, 7)
        
        flaimapper = FlaiMapper(args)
        
        i = 0
        for region in flaimapper.regions():
            if i == 0:
                self.assertEqual(region.region[0], 'chr1')
                
                self.assertEqual(region.region[1], 13-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 43+20+args.parameters.right_padding)
            elif i == 1:
                self.assertEqual(region.region[0], 'chr2')
                
                self.assertEqual(region.region[1], 54-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 54+26+args.parameters.right_padding)
            else:
                self.asserTrue(False,"Race condition?")
            
            i += 1
        
        self.assertEqual(i, 2)

    def test_01_c(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,'--verbose'])
        
        args.parameters.left_padding = 2
        args.parameters.right_padding = 2
        
        flaimapper = FlaiMapper(args)
        
        i = 0
        for region in flaimapper.regions():
            if i == 0:
                self.assertEqual(region.region[0], 'chr1')
                
                self.assertEqual(region.region[1], 13-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 13+20+4+args.parameters.right_padding)
            elif i == 1:
                self.assertEqual(region.region[0], 'chr1')
                
                self.assertEqual(region.region[1], 43-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 43+20+args.parameters.right_padding)
            elif i == 2:
                self.assertEqual(region.region[0], 'chr2')
                
                self.assertEqual(region.region[1], 54-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 54+26+args.parameters.right_padding)
            else:
                self.assertTrue(False,"Race condition?")
            i += 1
        
        self.assertEqual(i, 3)
    
    def test_01_d(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,'--verbose'])
        
        args.parameters.left_padding = 3
        args.parameters.right_padding = 3
        
        flaimapper = FlaiMapper(args)
        
        i = 0
        for region in flaimapper.regions():
            if i == 0:
                self.assertEqual(region.region[0], 'chr1')
                
                self.assertEqual(region.region[1], 13-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 13+20+4+args.parameters.right_padding)
            elif i == 1:
                self.assertEqual(region.region[0], 'chr1')
                
                self.assertEqual(region.region[1], 43-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 43+20+args.parameters.right_padding)
            else:
                self.assertEqual(region.region[0], 'chr2')
                
                self.assertEqual(region.region[1], 54-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 54+26+args.parameters.right_padding)
            
            i += 1
        
        self.assertEqual(i, 3)
    
    def test_01_e(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,'--verbose'])
        
        args.parameters.left_padding = 3
        args.parameters.right_padding = 4
        
        flaimapper = FlaiMapper(args)
        
        i = 0
        for region in flaimapper.regions():
            if i == 0:
                self.assertEqual(region.region[0], 'chr1')
                
                self.assertEqual(region.region[1], 13-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 43+20+args.parameters.right_padding)
            else:
                self.assertEqual(region.region[0], 'chr2')
                
                self.assertEqual(region.region[1], 54-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 54+26+args.parameters.right_padding)
            
            i += 1
        
        self.assertEqual(i, 2)
    
    def test_01_f(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,'--verbose'])
        
        args.parameters.left_padding = 4
        args.parameters.right_padding = 3
        
        flaimapper = FlaiMapper(args)
        
        i = 0
        for region in flaimapper.regions():
            if i == 0:
                self.assertEqual(region.region[0], 'chr1')
                
                self.assertEqual(region.region[1], 13-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 43+20+args.parameters.right_padding)
            else:
                self.assertEqual(region.region[0], 'chr2')
                
                self.assertEqual(region.region[1], 54-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 54+26+args.parameters.right_padding)
            
            i += 1
        
        self.assertEqual(i, 2)
    
    def test_01_g(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,'--verbose'])
        
        args.parameters.left_padding = 4
        args.parameters.right_padding = 4
        
        flaimapper = FlaiMapper(args)
        
        i = 0
        for region in flaimapper.regions():
            if i == 0:
                self.assertEqual(region.region[0], 'chr1')
                
                self.assertEqual(region.region[1], 13-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 43+20+args.parameters.right_padding)
            else:
                self.assertEqual(region.region[0], 'chr2')
                
                self.assertEqual(region.region[1], 54-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 54+26+args.parameters.right_padding)
            
            i += 1
        
        self.assertEqual(i, 2)
    
    
    def test_02(self):
        fname = 'test_FlaiMapper_test_02_output.gtf'
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,"-o",fname,'--verbose','--offset5p','4','--offset3p','4'])
        
        flaimapper = FlaiMapper(args)
        
        # Run analysis
        flaimapper.run()
        
        # assert Contents:
        self.assertTrue(filecmp.cmp(TESTS_FLAIMAPPER_TEST_02_OUTPUT_GTF, fname),msg="diff '"+TESTS_FLAIMAPPER_TEST_02_OUTPUT_GTF+"' '"+fname+"':\n"+subprocess.Popen(['diff',TESTS_FLAIMAPPER_TEST_02_OUTPUT_GTF,fname], stdout=subprocess.PIPE).stdout.read())
        
        os.remove(fname)# does not get executed if assertion is False

    def test_02b(self):
        fname = 'test_FlaiMapper_test_02b_output.gtf'
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,"-o",fname,'--verbose','--offset5p','5','--offset3p','5'])
        
        flaimapper = FlaiMapper(args)
        
        # Run analysis
        flaimapper.run()
        
        # assert Contents:
        assertion = (not filecmp.cmp(fname , TESTS_FLAIMAPPER_TEST_02_OUTPUT_GTF))
        self.assertTrue(assertion)
        
        if assertion:
            os.remove(fname)
    
    def test_03_a(self):
        fname = 'test_FlaiMapper_test_03_output.txt'
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,"-o",fname,"-f","1",'--verbose','--offset5p','4','--offset3p','4'])
        
        flaimapper = FlaiMapper(args)
        
        # Run analysis
        flaimapper.run()
        
        # assert Contents:
        self.assertTrue(filecmp.cmp(TESTS_FLAIMAPPER_TEST_03_a_OUTPUT_TXT, fname),msg="diff '"+TESTS_FLAIMAPPER_TEST_03_a_OUTPUT_TXT+"' '"+fname+"':\n"+subprocess.Popen(['diff',TESTS_FLAIMAPPER_TEST_03_a_OUTPUT_TXT,fname], stdout=subprocess.PIPE).stdout.read())

        os.remove(fname)
    
    def test_03_b(self):
        fname = 'test_FlaiMapper_test_03_fa_output.txt'
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,"-o",fname,"-f","1","--fasta",TESTS_FLAIMAPPER_FA,'--verbose','--offset5p','4','--offset3p','4'])
        
        flaimapper = FlaiMapper(args)
        
        # Run analysis
        flaimapper.run()
        
        # assert Contents:
        self.assertTrue(filecmp.cmp(TESTS_FLAIMAPPER_TEST_03_b_OUTPUT_TXT, fname),msg="diff '"+TESTS_FLAIMAPPER_TEST_03_b_OUTPUT_TXT+"' '"+fname+"':\n"+subprocess.Popen(['diff',TESTS_FLAIMAPPER_TEST_03_b_OUTPUT_TXT,fname], stdout=subprocess.PIPE).stdout.read())
        
        os.remove(fname)

    def test_04_a(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,'--verbose'])
        
        fm = FlaiMapper(args)
        
        # First check whether it's looking in the right regions
        i = 0
        for region in fm:
            if i == 0:
                self.assertEqual(region.region[0], 'chr1')
                self.assertEqual(region.region[1], max(0, 13-1-args.parameters.left_padding))
                self.assertEqual(region.region[2], 43+20+args.parameters.right_padding)
            elif i == 1:
                self.assertEqual(region.region[0], 'chr2')
                self.assertEqual(region.region[1], 54-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 54+26+args.parameters.right_padding)
            else:
                self.asserTrue(False,"Race condition?")
            i += 1
        self.assertEqual(i, 2)
        
        # Then detect the fragments
        i = 0
        for region in fm:
            for fragment in region:
                if i == 0:
                    self.assertEqual(region.region[1] + fragment.start, 13)
                    self.assertEqual(region.region[1] + fragment.stop, 36)
                elif i == 1:
                    self.assertEqual(region.region[1] + fragment.start, 43)
                    self.assertEqual(region.region[1] + fragment.stop, 62)
                elif i == 2:
                    self.assertEqual(region.region[1] + fragment.start, 54)
                    self.assertEqual(region.region[1] + fragment.stop, 79)
                else:
                    self.asserTrue(False,"Race condition?")
                i += 1
        self.assertEqual(i, 3)
    
    def test_04_b(self):
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,'-p',TESTS_FUNCTIONAL_DUCK7_PARAMS,'--verbose'])
        
        fm = FlaiMapper(args)
        
        # First check whether it's looking in the right regions
        i = 0
        for region in fm:
            if i == 0:
                self.assertEqual(region.region[0], 'chr1')
                
                self.assertEqual(region.region[1], 13-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 43+20+args.parameters.right_padding)
            elif i == 1:
                self.assertEqual(region.region[0], 'chr2')
                
                self.assertEqual(region.region[1], 54-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 54+26+args.parameters.right_padding)
            else:
                self.asserTrue(False,"Race condition?")
            i += 1
        
        self.assertEqual(i, 2)
        
        # Then detect the fragments
        i = 0
        for region in fm:
            for fragment in region:
                if i == 0:
                    self.assertEqual(region.region[1] + fragment.start, 13)
                    self.assertEqual(region.region[1] + fragment.stop, 36)
                elif i == 1:
                    self.assertEqual(region.region[1] + fragment.start, 43)
                    self.assertEqual(region.region[1] + fragment.stop, 62)
                elif i == 2:
                    self.assertEqual(region.region[1] + fragment.start, 54)
                    self.assertEqual(region.region[1] + fragment.stop, 79)
                else:
                    self.asserTrue(False,"Race condition?")
                i += 1
        
        self.assertEqual(i, 3)

    def test_04_c(self):
        """
        This should filter out the second fragment because it's 26bp away.
        """
        args = CLI([TESTS_EXAMPLE_ALIGNMENT_01,'-p',TESTS_FUNCTIONAL_DUCK26_PARAMS,'--verbose'])
        
        fm = FlaiMapper(args)
        
        # First check whether it's looking in the right regions
        i = 0
        for region in fm:
            if i == 0:
                self.assertEqual(region.region[0], 'chr1')
                
                self.assertEqual(region.region[1], max(0, 13-1-args.parameters.left_padding))
                self.assertEqual(region.region[2], 43+20+args.parameters.right_padding)
            elif i == 1:
                self.assertEqual(region.region[0], 'chr2')
                
                self.assertEqual(region.region[1], 54-1-args.parameters.left_padding)
                self.assertEqual(region.region[2], 54+26+args.parameters.right_padding)
            else:
                self.asserTrue(False,"Race condition?")
            i += 1
        self.assertEqual(i, 2)
        
        # Then detect the fragments
        i = 0
        for region in fm:
            for fragment in region:
                if i == 0:
                    self.assertEqual(region.region[1] + fragment.start, 13)
                    self.assertEqual(region.region[1] + fragment.stop, 36)
                elif i == 1:
                    self.assertEqual(region.region[1] + fragment.start, 54)
                    self.assertEqual(region.region[1] + fragment.stop, 79)
                else:
                    self.asserTrue(False,"Race condition?")
                i += 1
        
        self.assertEqual(i, 2)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
