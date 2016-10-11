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

import logging
from flaimapper.Data import *

class FilterParameters:
    """
    File format should be like this:

-3 \t 0.0
-2 \t 50.0
-1 \t 100.0
-1 \t 100.0
-2 \t 50.0
-3 \t 0.0

Every line represents an offset from a peak, and the second column
represents the percentage of reduction. The zero value is excluded.
    """
    def __init__(self,filename=None):
        self.parse(filename if filename != None else PARAMETERS_DEFAULT)
    
    def parse(self,filename):
        matrix = {}
        
        with open(filename, "r") as fh:
            for line in fh:
                line = line.strip()
                if len(line) > 0:
                    params = line.split("\t")
                    if len(params) != 2:
                        raise ValueError("inappropriate line detected in parameters file: %s:\n%s" % (filename, line))
                    
                    params = (int(params[0]),float(params[1]))
                    if params[0] != 0 and params[1] >= 0.0 and params[1] <= 100.0:
                        matrix[params[0]] = params[1]
                    else:
                        raise ValueError("inappropriate value(s) in line detected in parameters file %s" % (filename, line))
        
        self.set_matrix(matrix)
    
    def set_matrix(self,matrix):
        pos = []
        neg = []
        
        for key in matrix.keys():
            if key < 0:
                neg.append(key)
            else:
                pos.append(key)
        
        pos.sort()
        neg.sort()
        
        for i in range(1,len(pos)):
            if pos[i] - pos[i-1] != 1:
                raise ValueError("missing positive values in parameters file")

        for i in range(1,len(neg)):
            if neg[i] - neg[i-1] != 1:
                raise ValueError("missing negative values in parameters file")
        
        self.left_padding = len(pos)
        self.right_padding = len(neg)
        
        self.matrix = matrix
        logging.debug("Parsed filter parameters: [-"+str(self.left_padding)+","+str(self.right_padding)+"]")
