#!/usr/bin/env python

"""FlaiMapper: computational annotation of small ncRNA derived fragments using RNA-seq high throughput data

 Here we present Fragment Location Annotation Identification mapper
 (FlaiMapper), a method that extracts and annotates the locations of
 sncRNA-derived RNAs (sncdRNAs). These sncdRNAs are often detected in
 sequencing data and observed as fragments of their  precursor sncRNA.
 Using small RNA-seq read alignments, FlaiMapper is able to annotate
 fragments primarily by peak-detection on the start and  end position
 densities followed by filtering and a reconstruction processes.
 Copyright (C) 2011-2016:
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

from pkg_resources import resource_filename

PARAMETERS_DEFAULT = resource_filename("flaimapper","data/parameters.default.txt")

TESTS_EXAMPLE_ALIGNMENT_01 = resource_filename("flaimapper","data/tests/alignment_gap_6bp.bam")

TESTS_FLAIMAPPER_TEST_02_OUTPUT_GTF = resource_filename("flaimapper","data/tests/test_FlaiMapper_test_02_output.gtf")

TESTS_FLAIMAPPER_TEST_03_a_OUTPUT_TXT = resource_filename("flaimapper","data/tests/test_FlaiMapper_test_03_output.txt")
TESTS_FLAIMAPPER_TEST_03_b_OUTPUT_TXT = resource_filename("flaimapper","data/tests/test_FlaiMapper_test_03_fa_output.txt")

TESTS_FLAIMAPPER_FA = resource_filename("flaimapper","data/tests/test_FlaiMapper.fa")

TESTS_FUNCTIONAL_DUCK5_PARAMS = resource_filename("flaimapper","data/tests/test_functional.parameters.duck5.txt")
TESTS_FUNCTIONAL_DUCK6_PARAMS = resource_filename("flaimapper","data/tests/test_functional.parameters.duck6.txt")
TESTS_FUNCTIONAL_DUCK7_PARAMS = resource_filename("flaimapper","data/tests/test_functional.parameters.duck7.txt")
TESTS_FUNCTIONAL_DUCK15_PARAMS = resource_filename("flaimapper","data/tests/test_functional.parameters.duck15.txt")
TESTS_FUNCTIONAL_DUCK26_PARAMS = resource_filename("flaimapper","data/tests/test_functional.parameters.duck26.txt")
