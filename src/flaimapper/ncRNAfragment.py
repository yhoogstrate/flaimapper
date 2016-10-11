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

class ncRNAfragment:
    def __init__(self,start=None,stop=None,reference_sequence=None,masked_region=None,genomic_offset_masked_region=0):
        self.name = False
        
        self.set_location(start,stop,reference_sequence)
        self.genomic_offset_masked_region = genomic_offset_masked_region
        
        self.masked_region = masked_region
        self.sequence = None
        
        self.reset_supporting_reads()
        self.extended = {}
        self.evidence = None
    
    def get_sequence(self,fasta_handler=None):
        if(self.sequence):
            return self.sequence
        else:
            pass
    
    def set_sequence(self,sequence):
        self.sequence = sequence
    
    def set_name(self,name):
        self.name = name
    
    def set_evidence(self,evidence):
        self.evidence = evidence
    
    def get_name(self):
        return self.name
        
    def add_supporting_reads(self,number):
        self.supporting_reads += number
        
    def get_supporting_reads(self):
        return self.supporting_reads
    
    def spans_read(self,read,offset_left=0,offset_right=0):
        """
        no offset used..
        """
        return (((read.start+offset_left) >= self.start) and ((read.stop-offset_right) <= self.stop))
    
    def reset_supporting_reads(self):
        self.supporting_reads = 0										# all reads in-between the fragment
        self.supporting_reads_start = 0									# The reads with the start-position aligned exactly to the 5' of the fragment
        self.supporting_reads_stop = 0									# The reads with the end-position aligned exactly to the 3' of the fragment
    
    def set_location(self,start,stop,reference_sequence):
        self.start = start
        self.stop = stop
        self.reference_sequence = reference_sequence
    
    def get_start_position(self,absolute=False):
        if(absolute):
            self.start + self.genomic_offset_masked_region
        else:
            return self.start
    
    def get_stop_position(self,absolute=False):
        if(absolute):
            self.stop + self.genomic_offset_masked_region
        else:
            return self.stop
