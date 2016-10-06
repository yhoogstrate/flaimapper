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


import unittest,logging,os,subprocess
logging.basicConfig(format=flaimapper.__log_format__, level=logging.DEBUG)



class TestFunctional(unittest.TestCase):
    def test_01(self):
        pipe = subprocess.Popen(["flaimapper","-false-argument-"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout = pipe.stdout.read()
        stderr = pipe.stderr.read()
        exit_code = pipe.poll()

        self.assertTrue(exit_code != 0)
        self.assertTrue(stderr != '')
    
    def test_02(self):
        pipe = subprocess.Popen(["flaimapper","--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout = pipe.stdout.read()
        stderr = pipe.stderr.read()
        exit_code = pipe.poll()

        self.assertEqual(exit_code, 0)
        self.assertEqual(stderr, '')
    
    def test_03(self):
        # 01 extract BAM file:
        os.chdir("../share/small_RNA-seq_alignments/SRP028959/")
        subprocess.call(["tar","-xzf","SRR954958.tar.gz"])

def main():
    unittest.main()

if __name__ == '__main__':
    main()
