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
from flaimapper.utils import *


import unittest,logging,os,subprocess,shutil
logging.basicConfig(format=flaimapper.__log_format__, level=logging.DEBUG)



class TestFunctional(unittest.TestCase):
    def test_01(self):
        pipe = subprocess.Popen(["flaimapper","-false-argument-"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout = pipe.stdout.read()
        stderr = pipe.stderr.read()
        pipe.wait()
        exit_code = pipe.poll()

        self.assertTrue(exit_code != 0)
        self.assertTrue(stderr != '')
    
    def test_02(self):
        pipe = subprocess.Popen(["flaimapper","--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout = pipe.stdout.read()
        stderr = pipe.stderr.read()
        pipe.wait()
        exit_code = pipe.poll()

        self.assertEqual(exit_code, 0)
        self.assertEqual(stderr, '')
    
    def test_03(self):
        samplename = "SRR207111_HeLa18-30"
        
        os.chdir("../share/small_RNA-seq_alignments/SRP006788/")
        subprocess.call(["tar","-xzf",samplename+".tar.gz"])
        output_file = 'test.tabular.txt'
        
        pipe = subprocess.Popen(["flaimapper","-o",output_file,"-f","1",samplename+".bam"])
        pipe.wait()
        exit_code = pipe.poll()
        
        self.assertEqual(exit_code, 0)
        
        idx_test = parse_table(output_file)
        idx_old = parse_table('../../../output/FlaiMapper/SRP006788/01.a_output_flaimapper.txt',1,2,3,4)
        
        m = 0
        mm = 0
        mm_trna = 0
        """ finds new ones that were not found before"""
        for key in idx_test.keys():
            if key in idx_old.keys():
                if len(idx_test[key]) != len(idx_old[key]):
                    for key2 in idx_test[key].keys():
                        if key2 not in idx_old[key].keys():
                            if key.find('TRNA') == -1:
                                mm += 1
                            else:
                                mm_trna += 1
                else:
                    m += 1
            else:
                if key.find('TRNA') == -1:
                    mm += 1
                else:
                    mm_trna += 1
        
        total_evaluations = m+mm+mm_trna
        success = 1.0 - (1.0*(mm+mm_trna)/(total_evaluations))
        
        assertion1 = (mm <= 1)
        assertion2 = (mm_trna <= 63)
        assertion3 = (success >= 0.934)
        if assertion1 and assertion2 and assertion3:
            os.remove(output_file)
        
        os.remove(samplename+".bam")
        os.remove(samplename+".bam.bai")
        shutil.rmtree(samplename)
        
        self.assertTrue(assertion1,"Too many discrepancies with original results were found: %i non tRNA mispredictions" % mm)
        self.assertTrue(assertion2,"Too many discrepancies with original results were found: %i non tRNA mispredictions" % mm)# Latesr versions of FlaiMapper take into account the very last base as well, even if it's longer than the actual sequence length. This in particularly affected the results of tRNAs (due to CCA) 
        self.assertTrue(assertion3,"Too many discrepancies with original results were found; only: %d percent" % success)

    #def test_04
        #.... '-f','2'
        #parse_gff('test.gtf')

def main():
    unittest.main()

if __name__ == '__main__':
    main()
