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
import unittest
import logging
import os
import subprocess
import filecmp

from flaimapper.Data import TESTS_FUNCTIONAL_TEST_03_OUTPUT_TXT
from flaimapper.Data import TESTS_FUNCTIONAL_TEST_05
from flaimapper.Data import TESTS_FUNCTIONAL_TEST_05_OUTPUT_TXT
from flaimapper.utils import parse_table
from flaimapper.utils import parse_gff
from flaimapper.utils import get_file_diff

logging.basicConfig(format=flaimapper.__log_format__, level=logging.DEBUG)


class TestFunctional(unittest.TestCase):
    def test_01(self):
        with subprocess.Popen(["flaimapper", "-false-argument-"], stdout=subprocess.PIPE, stderr=subprocess.PIPE) as pipe:
            # stdout = pipe.stdout.read()
            stderr = pipe.stderr.read()

            pipe.wait()
            exit_code = pipe.poll()

            self.assertTrue(exit_code != 0)
            self.assertTrue(stderr != b'')

    def test_02(self):
        with subprocess.Popen(['flaimapper', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE) as pipe:
            # stdout = pipe.stdout.read()
            stderr = pipe.stderr.read()

            pipe.wait()
            exit_code = pipe.poll()
            pipe.kill()

            self.assertEqual(exit_code, 0)
            self.assertEqual(stderr, b'')

    def test_03(self):
        samplename = 'SRR207111_HeLa18-30'
        sampledir = '../share/small_RNA-seq_alignments/SRP006788/'

        command = ['tar', '-xzf', sampledir + samplename + '.tar.gz']
        subprocess.call(command, env=os.environ.copy())
        output_file = 'test_fuctional_03.tabular.txt'

        command = ['flaimapper', '-o', output_file, '-f', '1', '-r', '../share/annotations/ncRNA_annotation/ncrnadb09.fa', samplename + '.bam']
        with subprocess.Popen(command) as pipe:
            pipe.wait()
            exit_code = pipe.poll()

            self.assertEqual(exit_code, 0)

            idx_test = parse_table(output_file)
            idx_old = parse_table('../output/FlaiMapper/SRP006788/01.a_output_flaimapper.txt', 1, 2, 3, 4)

            m = 0
            mm = 0
            mm_trna = 0
            # finds new predictions that were not found in the results of the older versions
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
                        for key2 in idx_test[key]:
                            for item in idx_test[key][key2]:
                                self.assertTrue(item[3] != '', 'No sequence found for: %s' % item[0])  # No empty sequences
                        m += 1
                else:
                    if key.find('TRNA') == -1:
                        mm += 1
                    else:
                        mm_trna += 1

            total_evaluations = m + mm + mm_trna
            success = 1.0 - (1.0 * (mm + mm_trna) / total_evaluations)

            self.assertTrue(filecmp.cmp(TESTS_FUNCTIONAL_TEST_03_OUTPUT_TXT, output_file), msg="diff '" + TESTS_FUNCTIONAL_TEST_03_OUTPUT_TXT + "' '" + output_file + "':\n" + get_file_diff(TESTS_FUNCTIONAL_TEST_03_OUTPUT_TXT, output_file))

            assertion1 = (mm <= 1)
            assertion2 = (mm_trna <= 67)
            assertion3 = (success >= 0.9302)

            self.assertTrue(assertion1, "1. Too many discrepancies with original results were found: %i non tRNA mispredictions" % mm)
            self.assertTrue(assertion2, "2. Too many discrepancies with original results were found: %i tRNA mispredictions" % mm_trna)  # Latest versions of FlaiMapper take into account the very last base as well, even if it's longer than the actual sequence length. This in particularly affected the results of tRNAs (due to CCA)
            self.assertTrue(assertion3, "3. Too many discrepancies with original results were found; only: %f percent (m=%s, mm=%d, mm_trna=%d)" % (success, m, mm, mm_trna))

            # os.remove(samplename+".bam")
            # os.remove(samplename+".bam.bai")
            # shutil.rmtree(samplename)

    def test_04(self):
        samplename = 'SRR207111_HeLa18-30'
        sampledir = '../share/small_RNA-seq_alignments/SRP006788/'

        subprocess.call(['tar', '-xzf', sampledir + samplename + '.tar.gz'])
        output_file = 'test.gtf'

        command = ['flaimapper', '-o', output_file, '-f', '2', samplename + '.bam']
        with subprocess.Popen(command) as pipe:
            pipe.wait()
            exit_code = pipe.poll()

            self.assertEqual(exit_code, 0)
            data = parse_gff(output_file)
            u81_14 = False
            u81_46 = False
            u81_54 = False

            i = 0
            k = 0
            for chunk in data:  # range(len(data)):
                chunk = data[i]
                if chunk[0] == 'HGNC=10101&HUGO-Symbol=SNORD81&HUGO-Name=small_nucleolar_RNA,_C/D_box_81&LOCI=[chr1:173833274-173833370:strand=-]&SOURCE=RefSeq&SOURCE-ACCESSION=NR_003938&GENOME=hg19':
                    k += 1
                    if chunk[1] == 14 and chunk[2] == 36:
                        u81_14 = True
                    elif chunk[1] == 46 and chunk[2] == 67:
                        u81_46 = True
                    elif chunk[1] == 54 and chunk[2] == 79:
                        u81_54 = True
                i += 1

            if u81_14 and u81_46 and u81_54:
                os.remove(output_file)

            self.assertTrue(u81_14, "The first of the three SNORD81 fragments was not detected")
            self.assertTrue(u81_46, "The second of the three SNORD81 fragments was not detected")
            self.assertTrue(u81_54, "The third of the three SNORD81 fragments was not detected")
            self.assertTrue(k == (3 * 2), "More than 3 fragments (%i) of SNORD81 were detected" % k)  # *2 because every entry generates two GTF lines

            # os.remove(samplename+".bam")
            # os.remove(samplename+".bam.bai")
            # shutil.rmtree(samplename)

    def test_05(self):
        output_file = 'test.tabular.txt'

        command = ['flaimapper', '-o', output_file, '-f', '1', TESTS_FUNCTIONAL_TEST_05]
        with subprocess.Popen(command) as pipe:
            pipe.wait()
            exit_code = pipe.poll()

            self.assertEqual(exit_code, 0)

            # idx_test = parse_table(output_file)
            self.assertTrue(
                filecmp.cmp(TESTS_FUNCTIONAL_TEST_05_OUTPUT_TXT, output_file),
                msg="diff '" + TESTS_FUNCTIONAL_TEST_05_OUTPUT_TXT + "' '" + output_file + "':\n" + get_file_diff(TESTS_FUNCTIONAL_TEST_05_OUTPUT_TXT, output_file))

    def test_06(self):
        test_file_out = "tmp/test/test_functional_sslm2sam_test01.sam"
        test_file_cmp = "tests/results/test_functional_sslm2sam_test01.sam"

        command = ['sslm2sam', '-o', test_file_out, 'tests/data/sslm2sam']
        with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as pipe:
            # stdout = pipe.stdout.read()
            stderr = str(pipe.stderr.read())
            # pipe.wait()
            exit_code = pipe.returncode

            self.assertTrue(stderr != '', stderr)
            self.assertTrue(exit_code != 0, exit_code)

            self.assertTrue(filecmp.cmp(test_file_out, test_file_cmp), msg="diff '" + test_file_out + "' '" + test_file_cmp + "':\n" + get_file_diff(test_file_out, test_file_cmp))


def main():
    unittest.main()


if __name__ == '__main__':
    main()
