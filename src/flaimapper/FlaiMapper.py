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

import logging
import pysam

from .BAMParser import BAMParser
from .FragmentContainer import FragmentContainer
from .FragmentFinder import FragmentFinder


class FlaiMapper(FragmentContainer):
    def __init__(self,alignment_file):
        self.alignment = alignment_file
        self.sequences = {}
        
        logging.info(" - Initiated FlaiMapper Object")
        
        try:
            self.alignment.fetch()
        except:
            logging.info(' - Indexing BAM file with samtools: '+self.alignment.filename)
            pysam.index(self.alignment.filename)
            self.alignment = pysam.AlignmentFile(self.alignment.filename)
        
        try:
            self.alignment.fetch()
        except:
            raise Exception('Couldn\'t indexing BAM file with samtools: '+self.alignment.filename+'\nAre you sure samtools is installed?\n')
    
    def regions(self,filter_parameters):
        """
        Needs to find chunks of all consequently aligned blocks (+left
        and right padding distance of the filter)
        
        In case of large chromosome and filter of <-3,+3>:
        
        >chr1
        .....read.............................READ............
        .......read...........................................
        .........read..............................READ........
          [------------]                   [-------------]
        
        two regions to be yielded
        """
        
        i_dist_l = abs(filter_parameters.left_padding)
        i_dist_r = abs(filter_parameters.right_padding)
        i_dist = i_dist_l + i_dist_r
        
        for i in range(self.alignment.nreferences):
            s_name = self.alignment.references[i]
            ss = [None,None]
            
            for r in self.alignment.fetch(s_name):
                if len(r.blocks) > 0:
                    if ss[0] == None:
                        ss = [r.blocks[0][0],r.blocks[-1][1]-1]
                    else:
                        if r.blocks[-1][1]-1 > ss[1]:
                            if r.blocks[0][0] - ss[1] <= i_dist:
                                m  = max(ss[1],r.blocks[-1][1]-1)
                                ss[1] = m
                            else:
                                yield (s_name, max(0, ss[0] - i_dist_l - 1), max(0, ss[1] + i_dist_r + 1))
                                
                                ss = [r.blocks[0][0],r.blocks[-1][1]-1]
        
            if ss[0] != None:
                yield (s_name, max(0, ss[0] - i_dist_l - 1), max(0, ss[1] + i_dist_r + 1))
    
    def run(self,fasta_file,filter_parameters):
        logging.debug(" - Running fragment detection")
        self.fasta_file = fasta_file
        
        for region in self.regions(filter_parameters):
            logging.debug("   - Masked region: "+region[0]+":"+str(region[1])+"-"+str(region[2]))
            logging.debug("     * Acquiring statistics")
            
            # BAM
            aligned_reads = BAMParser(region,self.alignment)
            aligned_reads.parse_stats()
            logging.debug("     * Detecting fragments")
            
            fragments = FragmentFinder(aligned_reads, filter_parameters)
            fragments.run()
            self.add_fragments(fragments.results)
