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

import os,re,operator,argparse,sys,logging,subprocess

import pysam

from flaimapper.BAMParser import BAMParser
from flaimapper.SSLMParser import SSLMParser
from flaimapper.FragmentContainer import FragmentContainer
from flaimapper.FragmentFinder import FragmentFinder


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
            aligned_reads = BAMParser(region[0],region[1],region[2],self.alignment)
            aligned_reads.parse_stats()
            
            logging.debug("     * Detecting fragments")
            
            predicted_fragments = FragmentFinder(region, aligned_reads, filter_parameters, True)
            self.add_fragments(predicted_fragments, self.fasta_file)
    
    def count_reads_per_region(self,regions,links,masked_regions,reference_offset=0):
        """
        All sequences in our library of ncRNAs have been extended with 10 bases.
        """
        
        stats_table = {}
        stats_table['experimental']     = {'error_5p':{"<-5":0,-5:0,-4:0,-3:0,-2:0,-1:0,0:0,1:0,2:0,3:0,4:0,5:0,">5":0},'error_3p':{"<-5":0,-5:0,-4:0,-3:0,-2:0,-1:0,0:0,1:0,2:0,3:0,4:0,5:0,">5":0},'predicted':0,'not_predicted_no_reads':0,'not_predicted_with_reads':0}
        stats_table['not_experimental'] = {'error_5p':{"<-5":0,-5:0,-4:0,-3:0,-2:0,-1:0,0:0,1:0,2:0,3:0,4:0,5:0,">5":0},'error_3p':{"<-5":0,-5:0,-4:0,-3:0,-2:0,-1:0,0:0,1:0,2:0,3:0,4:0,5:0,">5":0},'predicted':0,'not_predicted_no_reads':0,'not_predicted_with_reads':0}
        
        logging.debug(" - Running fragment detection")
        
        i = 0
        j = 0
        
        #for ncRNA in self.alignment_directories_indexed.keys():
        for region in masked_regions:
            ncRNA = region[0]
            if(links.has_key(ncRNA)):
                logging.debug("   - Analysing: "+ncRNA)
                
                annotations = regions.index[links[ncRNA]]
                
                # Use SSLM2BAM if you need conversion - bam is the de facto standard
                aligned_reads = BAMParser(region[0],region[1],region[2],self.alignments)
                aligned_reads.parse_stats()
                
                predicted_fragments_obj = FragmentFinder(ncRNA,aligned_reads,None,True)
                predicted_fragments_obj.run()
                
                predicted_fragments = predicted_fragments_obj.results
                
                i += 1
                
                for annotation in annotations.fragments:
                    closest = self.find_closest_overlapping_fragment(annotation,predicted_fragments,reference_offset)
                    j += 1
                    
                    if(closest):
                        errors = self.find_errors(annotation,closest,reference_offset)
                        err_5p = errors[0]
                        err_3p = errors[1]
                        
                        if(err_5p > 5):
                            err_5p = ">5"
                        elif(err_5p < -5):
                            err_5p = "<-5"
                        
                        if(err_3p > 5):
                            err_3p = ">5"
                        elif(err_3p < -5):
                            err_3p = "<-5"
                        
                        if(annotation.evidence == "experimental"):
                            stats_table['experimental']["predicted"] += 1
                            stats_table['experimental']["error_5p"][err_5p] += 1
                            stats_table['experimental']["error_3p"][err_3p] += 1
                            
                        else:
                            stats_table['not_experimental']["predicted"] += 1
                            stats_table['not_experimental']["error_5p"][err_5p] += 1
                            stats_table['not_experimental']["error_3p"][err_3p] += 1
                        
                    else:
                        if(annotation.evidence == "experimental"):
                            if(annotation.get_supporting_reads() == 0):
                                stats_table['experimental']["not_predicted_no_reads"] += 1
                            else:
                                stats_table['experimental']["not_predicted_with_reads"] += 1
                        else:
                            if(annotation.get_supporting_reads() == 0):
                                stats_table['not_experimental']["not_predicted_no_reads"] += 1
                            else:
                                stats_table['not_experimental']["not_predicted_with_reads"] += 1
        
        logging.info(str(i)+" annotated pre-miRNAs")
        logging.info(str(j)+" annotated miRNAs")
        
        return stats_table
    
    def count_error_with_intensity(self,regions,links,masked_regions,reference_offset=0):
        """All sequences in our library of ncRNAs were extended with 10 bases.
        """
        out = []
        logging.debug(" - Running fragment detection")
        
        for region in masked_regions:
            ncRNA = region[0]
            if(links.has_key(ncRNA)):
                logging.debug("   - Analysing: "+ncRNA)
                
                annotations = regions.index[links[ncRNA]]
                
                # Bam only
                aligned_reads = BAMParser(region[0],region[1],region[2],self.alignments)
                aligned_reads.parse_stats()
                
                predicted_fragments_obj = FragmentFinder(ncRNA,aligned_reads)
                predicted_fragments_obj.run()
                predicted_fragments = predicted_fragments_obj.results
                
                aligned_reads.count_reads_per_region(predicted_fragments_obj.results)
                
                for mirna_annotation in annotations.fragments:
                    closest_fragment = self.find_closest_overlapping_fragment(mirna_annotation,predicted_fragments,reference_offset)
                    
                    if(closest_fragment):
                        errors = self.find_errors(mirna_annotation,[(closest_fragment.start - reference_offset), (closest_fragment.stop - reference_offset)])#@todo ,reference_offset
                        err_5p = errors[0]
                        err_3p = errors[1]
                        
                        #out.append({'5p':[closest_fragment[2],err_5p],'3p':[closest_fragment[3],err_3p]})
                        out.append({'5p':[closest_fragment.supporting_reads_start,err_5p],'3p':[closest_fragment.supporting_reads_stop,err_3p],'coverage':closest_fragment.supporting_reads})
        
        return out
    
    def find_closest_overlapping_fragment(self,annotated_fragment,predicted_fragments,reference_offset=0):
        closest = False
        closest_overlapping_bases = 0
        for predicted_fragment in predicted_fragments:
            #predicted = [(predicted_fragment["start"] - reference_offset), (predicted_fragment["stop"] - reference_offset),predicted_fragment.start_supporting_reads,predicted_fragment.stop_supporting_reads]
            overlap = self.find_overlapping_bases([annotated_fragment.start,annotated_fragment.stop],[(predicted_fragment["start"] - reference_offset), (predicted_fragment["stop"] - reference_offset)])
            if(overlap > 0 and overlap > closest_overlapping_bases):
                closest_overlapping_bases = overlap
                closest = predicted_fragment
        
        return closest
    
    def find_overlapping_bases(self,fragment_1,fragment_2):
        if(fragment_2[0] < fragment_1[0]):
            return self.find_overlapping_bases(fragment_2,fragment_1)
        else:
            return fragment_1[1] - fragment_2[0]
    
    def find_errors(self,annotated_fragment,predicted_fragment,reference_offset=0):
        """
        Example:
               [ miRNA ]
            [ fragment* ]
        
        mirna: 4,12
        fragment: 0,14
        
        error_5p = 0 - 4 = -4
        error_3p = 13 - 12 = 1
        
        """
        
        error_5p = predicted_fragment.start - annotated_fragment.start - reference_offset
        error_3p = predicted_fragment.stop - annotated_fragment.stop - reference_offset
        
        return [error_5p,error_3p]
