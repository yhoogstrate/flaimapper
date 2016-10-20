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

from flaimapper.FragmentFinder import FragmentFinder
from flaimapper.FilterParameters import FilterParameters
from flaimapper.MaskedRegion import MaskedRegion
from flaimapper.Data import *
import flaimapper
import unittest,logging
logging.basicConfig(format=flaimapper.__log_format__, level=logging.DEBUG)


class TestFragmentFinder(unittest.TestCase):
    def test_01(self):
        """Test alignment:
        
        ..............
         [----]
         [----]
               [----]
        """
        matrices = MaskedRegion(('chr1',0,100))

        matrices.start_positions =   [0,2,0,0,0,0, 0,1,0,0,0,0, 0,0]
        matrices.stop_positions =    [0,0,0,0,0,0, 2,0,0,0,0,0, 1,0]
        
        matrices.start_avg_lengths = [0,6,0,0,0,0, 0,6,0,0,0,0, 0,0]
        matrices.stop_avg_lengths =  [0,0,0,0,0,0,-6,0,0,0,0,0,-6,0]
        
        # for the secoond fragment duck5 should work, duck6 should not.
        fp5 = FilterParameters(TESTS_FUNCTIONAL_DUCK5_PARAMS)
        ff5 = FragmentFinder(matrices,fp5)
        ff5.run()
        
        self.assertEqual(len(ff5.results), 2)
        
        self.assertEqual(ff5.results[0].start, 1)
        self.assertEqual(ff5.results[0].stop, 6)

        self.assertEqual(ff5.results[1].start, 7)
        self.assertEqual(ff5.results[1].stop, 12)
        
        
        fp6 = FilterParameters(TESTS_FUNCTIONAL_DUCK6_PARAMS)
        ff6 = FragmentFinder(matrices,fp6)
        ff6.run()

        self.assertEqual(len(ff6.results), 1)
        
        self.assertEqual(ff6.results[0].start, 1)
        self.assertEqual(ff6.results[0].stop, 6)

        
        fp7 = FilterParameters(TESTS_FUNCTIONAL_DUCK7_PARAMS)
        ff7 = FragmentFinder(matrices,fp7)
        ff7.run()
        
        self.assertEqual(len(ff7.results), 1)
        
        self.assertEqual(ff7.results[0].start, 1)
        self.assertEqual(ff7.results[0].stop, 6)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
