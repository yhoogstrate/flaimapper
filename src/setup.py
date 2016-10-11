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

from distutils.core import setup
from setuptools import setup, find_packages

setup(name='flaimapper',
        version=flaimapper.__version__,
        description='Fragment Location Annotation Identification Mapper',
        author=flaimapper.__author__,
        maintainer=flaimapper.__author__,
        url='https://github.com/yhoogstrate/flaimapper',
        
        scripts=["bin/flaimapper","bin/sslm2sam"],
        packages=['flaimapper'],
        #package_dir={'flaimapper': 'flaimapper'},
        package_data={'': ['data/*.*','data/tests/*.*']},
        include_package_data=True,
        
        # Very severe backwards incompatibility in 0.9 and above
        setup_requires=['pysam >= 0.8.0,<0.9','nose2'],
        install_requires=['pysam >= 0.8.0,<0.9'],
        
        test_suite="tests",
        
        classifiers=[
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Operating System :: OS Independent',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            ],
    )
