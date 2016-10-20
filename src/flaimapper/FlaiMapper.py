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

import flaimapper
from .BAMParser import BAMParser
from .FragmentFinder import FragmentFinder


class FlaiMapper():
    def __init__(self,settings):
        self.settings = settings
        self.settings.alignment_file = self.settings.alignment_file
        self.sequences = {}
        
        logging.info(" - Initiated FlaiMapper Object")
        
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

    def __iter__(self):
        for uid in sorted(self.sequences.keys()):
            for reference_sequence in self.sequences[uid]:
                for fragment in reference_sequence:
                    yield fragment
    
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
                                yield (s_name, max(0, ss[0] - i_dist_l - 1), max(0, ss[1] + i_dist_r + 1))
                                
                                ss = [r.blocks[0][0],r.blocks[-1][1]-1]
        
            if ss[0] != None:
                yield (s_name, max(0, ss[0] - i_dist_l - 1), max(0, ss[1] + i_dist_r + 1))
    
    def run(self):#,fasta_file,filter_parameters):
        logging.debug(" - Running fragment detection")
        
        i = 0
        for region in self.regions(self.settings.parameters):#@todo remove argument - it's part of self/settings anyway
            logging.debug("   - Masked region: "+region[0]+":"+str(region[1])+"-"+str(region[2]))
            logging.debug("     * Acquiring statistics")
            
            # BAM
            aligned_reads = BAMParser(region,self.settings.alignment_file)
            aligned_reads.parse_stats()
            logging.debug("     * Detecting fragments")
            
            fragments = FragmentFinder(aligned_reads, self.settings.parameters)
            fragments.run()
            self.add_fragments(fragments.results)
            i += len(fragments.results)
        
        logging.info(' - Detected %i fragments' % i)

    
    #@tofo get rid of inserting fasta_file HERE 
    def add_fragments(self,fragment_finder_results):
        #for fragment in fragment_finder_results.results:
        if len(fragment_finder_results) > 0:
            region = fragment_finder_results[0].masked_region
            uid = region[0]+"_"+str(region[1])+"_"+str(region[2])
            if(uid not in self.sequences.keys()):
                self.sequences[uid] = []
            
            self.sequences[uid].append(fragment_finder_results)
    
    def export_table(self,filename):
        """Exports the discovered fragments to a tab-delimited file.
        
        The following format is exported:
        
        ----
        @param filename The target file.
        
        @return:
        @rtype:
        """
        if(not self.sequences):
            logging.warning("     * Warning: no fragments detected")
        else:
            if(filename == "-"):
                fh = sys.stdout
            else:
                fh = open(filename,'w')
            
            if(self.settings.fasta_handle):
                fh.write("Fragment\tSize\tReference sequence\tStart\tEnd\tPrecursor\tStart in precursor\tEnd in precursor\tSequence (no fasta file given)\tCorresponding-reads (start)\tCorresponding-reads (end)\tCorresponding-reads (total)\n")
            else:
                fh.write("Fragment\tSize\tReference sequence\tStart\tEnd\tPrecursor\tStart in precursor\tEnd in precursor\tSequence\tCorresponding-reads (start)\tCorresponding-reads (end)\tCorresponding-reads (total)\n")
            
            for uid in sorted(self.sequences.keys()):
                for result in self.sequences[uid]:
                    if result:
                        name = result[0].masked_region[0]
                        fragments_sorted_keys = {}
                        for fragment in result:
                            fragments_sorted_keys[fragment.start] = fragment
                        
                        i = 0
                        for key in sorted(fragments_sorted_keys.keys()):	# Walk over i in the for-loop:
                            i += 1
                            fragment = fragments_sorted_keys[key]
                            fragment_uid = 'FM_'+result[0].masked_region[0]+'_'+str(i).zfill(12)
                            fh.write(fragment.to_table_entry(fragment_uid, result[0].masked_region, self.settings.fasta_handle))
            
            fh.close()
    
    def open_gtf(self,filename):
        if(filename == "-"):
            fh = sys.stdout
        else:
            fh = open(filename,'w')
        
        return fh
    
    def export_gtf(self,filename,offset5p,offset3p):
        if(filename == "-"):
            fh = sys.stdout
        else:
            fh = open(filename,'w')
        
        for uid in sorted(self.sequences.keys()):
            for result in self.sequences[uid]:
                if result:
                    name = result[0].masked_region[0]
                    fragments_sorted_keys = {}
                    for fragment in result:
                        fragments_sorted_keys[fragment.start] = fragment
                    
                    i = 0
                    for key in sorted(fragments_sorted_keys.keys()):# Walk over i in the for-loop:
                        i += 1
                        fragment = fragments_sorted_keys[key]
                        fragment_uid = 'FM_'+result[0].masked_region[0]+'_'+str(i).zfill(12)
                        fh.write(fragment.to_gtf_entry(fragment_uid, result[0].masked_region, offset5p, offset3p))
        
        fh.close()
    
    def write(self,export_format,output_filename,offset5p,offset3p):
        logging.debug(" - Exporting results to: "+output_filename)
        
        if(export_format == 1):
            logging.info("   - Format: tab-delimited, per fragment")
            self.export_table(output_filename)
        elif(export_format == 2):
            logging.info("   - Format: GTF")
            self.export_gtf(output_filename,offset5p,offset3p)
