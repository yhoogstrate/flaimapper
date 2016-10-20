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
import flaimapper


class FragmentContainer():
    def __iter__(self):
        for uid in sorted(self.sequences.keys()):
            for reference_sequence in self.sequences[uid]:
                for fragment in reference_sequence:
                    yield fragment
    
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
            
            
            if(self.fasta_file):
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
                            
                            # Fragment uid
                            fh.write('FM_'+ name+'_'+str(i).zfill(12)+"\t")
                            
                            # Size
                            fh.write(str(fragment.stop - fragment.start + 1) + "\t")
                            
                            # Reference sequence 
                            fh.write(name + "\t")
                            
                            # Start
                            fh.write(str(fragment.start) + "\t")
                            
                            # End
                            fh.write(str(fragment.stop)+"\t")
                            
                            # Precursor
                            fh.write(name)
                            
                            # Start in precursor
                            fh.write("\t" + str(fragment.start-fragment.masked_region[1])+ "\t")
                            
                            # End in precursor
                            fh.write(str(fragment.stop-fragment.masked_region[1])+"\t")
                            
                            # Sequence 
                            if(self.fasta_file):
                                # PySam 0.8.2 claims to use 0-based coordinates pysam.FastaFile.fetch().
                                # This is only true for the start position, the end-position is 1-based.
                                fh.write(str(self.fasta_file.fetch(name,fragment.start,fragment.stop+1)))
                            
                            # Start supporting reads
                            fh.write("\t"+str(fragment.supporting_reads_start)+"\t")
                            
                            # Stop supporting reads
                            fh.write(str(fragment.supporting_reads_stop)+"\t")
                            
                            # Total supporting reads
                            fh.write(str(fragment.supporting_reads_stop+fragment.supporting_reads_start) + "\n")
            
            fh.close()
    
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
                        
                        ## Line 1: type sncdRNA
                        # Seq-name
                        fh.write(name + "\t")
                        
                        # Source
                        fh.write("flaimapper-v"+flaimapper.__version__+"\t")
                        
                        # Feature
                        fh.write("sncdRNA\t")
                        
                        # Start
                        fh.write(str(fragment.start+1) + "\t")
                        
                        # End
                        fh.write(str(fragment.stop+1)+"\t")
                        
                        # Score
                        fh.write(str(fragment.supporting_reads_stop+fragment.supporting_reads_start) + "\t")
                        
                        # Strand and Frame
                        fh.write(".\t.\t")
                        
                        # Attribute
                        attributes = []
                        attributes.append('gene_id "FM_'+ name+'_'+str(i).zfill(12)+'"' )
                        
                        fh.write(", ".join(attributes)+"\n")
                        
                        
                        ## Line 2: type exon, with offset used for counting in e.g. HTSeq-count / featureCounts
                        
                        ## Line 1: type sncdRNA
                        # Seq-name
                        fh.write(name + "\t")
                        
                        # Source
                        fh.write("flaimapper-v"+flaimapper.__version__+"\t")
                        
                        # Feature
                        fh.write("exon\t")
                        
                        # Start
                        fh.write(str(max(1,fragment.start+1-offset5p))+"\t")
                        
                        # End
                        fh.write(str(max(1,fragment.stop+1+offset3p))+"\t")
                        
                        # Score
                        fh.write(str(fragment.supporting_reads_stop+fragment.supporting_reads_start) + "\t")
                        
                        # Strand and Frame
                        fh.write(".\t.\t")
                        
                        # Attribute
                        attributes = []
                        attributes.append('gene_id "FM_'+ name+'_'+str(i).zfill(12)+'"' )
                        
                        fh.write(", ".join(attributes)+"\n")
                        
        fh.close()
    
    def write(self,export_format,output_filename,offset5p,offset3p):
        logging.debug(" - Exporting results to: "+output_filename)
        
        if(export_format == 1):
            logging.info("   - Format: tab-delimited, per fragment")
            self.export_table(output_filename)
        elif(export_format == 2):
            logging.info("   - Format: GTF")
            self.export_gtf(output_filename,offset5p,offset3p)
