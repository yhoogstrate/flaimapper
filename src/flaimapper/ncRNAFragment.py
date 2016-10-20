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


class ncRNAFragment:
    def __init__(self,reference,start,stop):
        self.reference = reference
        self.start = start
        self.stop = stop
        
        self.supporting_reads = 0			# all reads in-between the fragment
        self.supporting_reads_start = 0		# The reads with the start-position aligned exactly to the 5' of the fragment
        self.supporting_reads_stop = 0		# The reads with the end-position aligned exactly to the 3' of the fragment
    
    def to_gtf_entry(self, uid, type_exon_offset5p, type_exon_offset3p):
        ## Line 1: type sncdRNA
        out_str = ("%s\t" # Reference
                   "flaimapper-v%s\t" # Source
                   "sncdRNA\t"
                   "%i\t"# Start
                   "%i\t"# End
                   "%i\t"# Score
                   ".\t.\t"# Strand and Frame
                   'gene_id "%s"\n' # Attributes (gene_id only)
                   ) % (self.reference,
                        flaimapper.__version__,
                        self.start+1,
                        self.stop+1,
                        self.supporting_reads_stop+self.supporting_reads_start,
                        uid)
        
        ## Line 2: type exon, with offset used for counting in e.g. HTSeq-count / featureCounts
        out_str += ("%s\t" # Reference
                    "flaimapper-v%s\t" # Source
                    "exon\t"
                    "%i\t"# Start
                    "%i\t"# End
                    "%i\t"# Score
                    ".\t.\t"# Strand and Frame
                    'gene_id "%s"\n' # Attributes (gene_id only)
                    ) % (self.reference,
                         flaimapper.__version__,
                         max(1,self.start+1-type_exon_offset5p),
                         max(1,self.stop+1+type_exon_offset3p),
                         self.supporting_reads_stop+self.supporting_reads_start,
                         uid)
        
        return out_str
    
    def to_table_entry(self, uid, masked_region, fasta_file):
        return ("%s\t"
                   "%i\t"
                   "%s\t"
                   "%i\t"
                   "%i\t"
                   "%s\t"
                   "%i\t"
                   "%i\t"
                   "%s\t"
                   "%i\t"
                   "%i\t"
                   "%i\n"
                   ) % (uid,
                        self.stop - self.start + 1,
                        self.reference,
                        self.start,
                        self.stop,
                        self.reference,
                        self.start-masked_region.region[1],
                        self.stop-masked_region.region[1],
                        fasta_file.fetch(self.reference,self.start,self.stop+1) if fasta_file else '',
                        self.supporting_reads_start,
                        self.supporting_reads_stop,
                        self.supporting_reads_stop+self.supporting_reads_start)
