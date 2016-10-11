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

from flaimapper.FilterParameters import FilterParameters
from flaimapper.Data import *
import unittest,logging
logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)


class TestFilterParameters(unittest.TestCase):
    def test_01(self):
        fp = FilterParameters()
        
        self.assertEqual(fp.left_padding, 15)
        self.assertEqual(fp.right_padding, 15)

        keys = fp.matrix.keys()
        for i in range(-15,0):
            self.assertTrue(i in keys)
            self.assertTrue(fp.matrix[i] >= 0.0)
            self.assertTrue(fp.matrix[i] <= 100.0)
            
    def test_02(self):
        fp = FilterParameters(TESTS_FUNCTIONAL_DUCK7_PARAMS)
        
        self.assertEqual(fp.left_padding, 7)
        self.assertEqual(fp.right_padding, 7)
        
        keys = fp.matrix.keys()
        for i in range(-7,0):
            self.assertTrue(i in keys)
            self.assertTrue(fp.matrix[i] >= 0.0)
            self.assertTrue(fp.matrix[i] <= 100.0)

def main():
    unittest.main()

if __name__ == '__main__':
    main()
