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
import unittest,subprocess,os,shutil,logging

logging.basicConfig(format=flaimapper.__log_format__, level=logging.DEBUG)

from flaimapper.CLI import CLI_sslm2sam
from flaimapper.SSLMParser import SSLMParser


class TestFlaiMapper(unittest.TestCase):
    def test_01_a(self):
        subprocess.call(["tar","-xzf","../share/small_RNA-seq_alignments/SRP028959/SRR954958.tar.gz"])
        
        args = CLI_sslm2sam(['-o','test.sam','SRR954958'])
        sslm2bed_converter = SSLMParser(args.sslm_directory)
        sslm2bed_converter.convert_to_sam(args.output)
        
        assertion = (os.stat("test.sam").st_size == 46985661)
        self.assertTrue(assertion, "Incorrect ../share/small_RNA-seq_alignments/SRP028959/test.sam")# Assume file size is sufficient :)
        
        if assertion:
            os.remove("test.sam")
        
        os.remove("SRR954958.bam")
        os.remove("SRR954958.bam.bai")
        shutil.rmtree("SRR954958")


def main():
    unittest.main()

# Tests fall outside __main__
if __name__ == '__main__':
    main()
