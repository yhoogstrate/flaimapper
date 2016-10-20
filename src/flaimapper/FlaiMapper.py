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
import logging,sys
logging.basicConfig(format=flaimapper.__log_format__, level=logging.DEBUG)

import pysam

from .MaskedRegion import MaskedRegion


class FlaiMapper():
    def __init__(self,settings):
        self.settings = settings
        
        logging.info(" - Initiated FlaiMapper Object")
        self.check_alignment_index()
    
    def check_alignment_index(self):
        try:
            self.settings.alignment_file.fetch()
        except:
            logging.info(' - Indexing BAM file with samtools: '+self.settings.alignment_file.filename)
            pysam.index(self.settings.alignment_file.filename)
            self.settings.alignment_file = pysam.AlignmentFile(self.settings.alignment_file.filename)
        
        try:
            self.settings.alignment_file.fetch()
        except:
            raise Exception('Couldn\'t indexing BAM file with samtools: '+self.settings.alignment_file.filename+'\nAre you sure samtools is installed?\n')
    
    def regions(self):
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
        
        i_dist_l = abs(self.settings.parameters.left_padding)
        i_dist_r = abs(self.settings.parameters.right_padding)
        i_dist = i_dist_l + i_dist_r
        
        for i in range(self.settings.alignment_file.nreferences):
            s_name = self.settings.alignment_file.references[i]
            ss = [None,None]
            
            for r in self.settings.alignment_file.fetch(s_name):
                if len(r.blocks) > 0:
                    if ss[0] == None:
                        ss = [r.blocks[0][0],r.blocks[-1][1]-1]
                    else:
                        if r.blocks[-1][1]-1 > ss[1]:
                            if r.blocks[0][0] - ss[1] <= i_dist:
                                m  = max(ss[1],r.blocks[-1][1]-1)
                                ss[1] = m
                            else:
                                yield MaskedRegion((s_name, max(0, ss[0] - i_dist_l - 1), max(0, ss[1] + i_dist_r + 1)),self.settings)
                                
                                ss = [r.blocks[0][0],r.blocks[-1][1]-1]
            
            if ss[0] != None:
                yield MaskedRegion((s_name, max(0, ss[0] - i_dist_l - 1), max(0, ss[1] + i_dist_r + 1)),self.settings)
    
    def __iter__(self):
        for region in self.regions():
            yield region
    
    def run(self):
        if(self.settings.format == 1):
            fh = self.open_table()
        elif(self.settings.format == 2):
            fh = self.open_gtf()
        
        logging.debug(" - Starting fragment detection")
        
        k = 0
        for region in self:
            i = 0
            for fragment in region:
                i += 1
                fragment_uid = 'FM_'+region.region[0]+'_'+str(i).zfill(12)
                
                if(self.settings.format == 1):
                    fh.write(fragment.to_table_entry(fragment_uid, region, self.settings.fasta_handle))
                elif(self.settings.format == 2):
                    fh.write(fragment.to_gtf_entry(fragment_uid, self.settings.offset5p, self.settings.offset3p))
            
            k += i
        
        fh.close()
        logging.info(' - Detected %i fragments' % k)

    def open_gtf(self):
        logging.info(" - Exporting results to: "+self.settings.output+" (GTF)")
        
        if(self.settings.output == "-"):
            fh = sys.stdout
        else:
            fh = open(self.settings.output,'w')
        
        return fh
    
    def open_table(self):
        logging.info(" - Exporting results to: "+self.settings.output+" (tab-delimited, per fragment)")
        
        if(self.settings.output == "-"):
            fh = sys.stdout
        else:
            fh = open(self.settings.output,'w')
        
        if(self.settings.fasta_handle):
            fh.write("Fragment\tSize\tReference sequence\tStart\tEnd\tPrecursor\tStart in precursor\tEnd in precursor\tSequence (no fasta file given)\tCorresponding-reads (start)\tCorresponding-reads (end)\tCorresponding-reads (total)\n")
        else:
            fh.write("Fragment\tSize\tReference sequence\tStart\tEnd\tPrecursor\tStart in precursor\tEnd in precursor\tSequence\tCorresponding-reads (start)\tCorresponding-reads (end)\tCorresponding-reads (total)\n")
        
        return fh
