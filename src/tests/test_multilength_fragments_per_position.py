#!/usr/bin/env python

"""FlaiMapper: computational annotation of small ncRNA derived fragments using RNA-seq high throughput data

 Here we present Fragment Location Annotation Identification mapper
 (FlaiMapper), a method that extracts and annotates the locations of
 sncRNA-derived RNAs (sncdRNAs). These sncdRNAs are often detected in
 sequencing data and observed as fragments of their  precursor sncRNA.
 Using small RNA-seq read alignments, FlaiMapper is able to annotate
 fragments primarily by peak-detection on the start and  end position
 densities followed by filtering and a reconstruction processes.
 Copyright (C) 2011-2018:
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
import unittest
import os
import logging
import pysam

from flaimapper.FlaiMapper import FlaiMapper
from flaimapper.CLI import CLI


logging.basicConfig(format=flaimapper.__log_format__, level=logging.DEBUG)


class TestFlaiMapper3(unittest.TestCase):
    def test_01(self):
        basename = 'multilength_fragments_per_position_001'

        if not os.path.exists('tmp/' + basename + '.bam'):
            fhq = open('tmp/' + basename + '.bam', "wb")
            fhq.write(pysam.view('-bS', 'tests/data/' + basename + ".sam"))
            fhq.close()

        args = CLI(['tmp/' + basename + '.bam', '--verbose'])

        args.parameters.left_padding = 0
        args.parameters.right_padding = 0

        flaimapper = FlaiMapper(args)
        i = 0
        for region in flaimapper.regions():
            self.assertEqual(region.region[0], 'SNORD78')
            for result in region:
                if i == 0:
                    self.assertEqual(region.region[1] + result.start, 11)
                    self.assertEqual(region.region[1] + result.stop, 11 + 61)
                elif i == 1:
                    self.assertEqual(region.region[1] + result.start, 44)
                    self.assertEqual(region.region[1] + result.stop, 44 + 28)
                i += 1
        self.assertEqual(i, 2)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
