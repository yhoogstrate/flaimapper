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


import os,re,operator,argparse,sys,logging


from flaimapper.Read import Read
from flaimapper.ncRNAfragment import ncRNAfragment
from flaimapper.MaskedRegion import MaskedRegion


class SSLMParser():
    """parseNcRNA is a class that parses the SSLM alignment files.
    """
    regex1 = re.compile("^>(.*?)_x([0-9]+)$")
    
    def __init__(self, sslm_directory):
        self.sslm_directory = sslm_directory
    
    def get_length(self,filename):
        i = 0
        with open(filename,"rU") as fh:
            for line in fh:
                if i == 1:
                    line = line.strip()
                    return len(line)
                i += 1
        
        raise Exception("File "+filename+" did not contain a seuqence that can be used to estimate length")
    
    def parse_reads(self,filename):
        """parse the reads from a SSLM (FASTA) file and return each read
        as an iterator object
        """
        
        previous_line = ""
        
        i = 0
        with open(filename,'r') as fh:
            for line in fh:
                line = line.strip().replace('revcomp','')
                
                if(i % 2 == 1):
                    if(i == 1):
                        self.sequence = line
                    else:
                        # previous_line = ">fasta name _hits etc"
                        #          line = "-----ACTG-----"
                        
                        k = previous_line.lower().find('_hits')
                        
                        if(k > -1):
                            name = previous_line[1:k]
                            numberofhits = int(previous_line[k+5::])
                        else:
                            m = self.regex1.search(name)			# For the "_x123" suffix
                            
                            if(m):
                                name = m.group(1)
                                numberofhits = int(m.group(2))
                            else:
                                name = previous_line[::-1].lstrip(">")
                                numberofhits = 1
                        
                        start_pos = self.get_start_position(line)
                        stop_pos = self.get_stop_position(line)
                        
                        for j in range(numberofhits):
                            yield Read(start_pos,stop_pos,name,line[start_pos:stop_pos])
                else:
                    previous_line = line
                
                i += 1
    
    def get_alignment_files(self):
        idx = {}
        with open(self.sslm_directory+"/idreadable.txt",'rU') as fh:
            for line in fh:
                line = line.strip()
                if line not in ["","sequence\tfilename"]:
                    line = line.split("\t")
                    idx[line[0][1:]] = self.sslm_directory+"/validated/"+line[1]+".fa"
        
        for key in sorted(idx.keys()):
            yield key,idx[key]
    
    def get_start_position(self,read,extention='-'):
        """
        Finds start-position of a read according to lines of the following format:
        -----TACCCTGTAGAGCCGAATTTGT-----
             *
        
        This example will return: 5 because at the 5th position is a 'T' (notice that we count from 0)
        
        ----
        @return:
        @rtype: integer
        """
        
        return len(read)-len(read.lstrip(extention))
    
    def get_stop_position(self,read,extention='-'):
        """
        Finds stop-position of a read according to lines of the following format:
        -----TACCCTGTAGAGCCGAATTTGT-----
                                   *
        
        This example will return: 27 because at the 27th position is a '-' (notice that we count from 0)
        
        ----
        @return:
        @rtype: integer
        """
        
        return len(read.rstrip(extention))
    
    def parse_regions(self):
        for region in self.get_alignment_files():
            yield region[0],0,self.get_length(region[1])-1,region[1]# zero based: if len == 1, coord will be [0, 0]
    
    def convert_to_sam(self,output):
        logging.debug("   - Converting to SAM: "+output)
        
        if(output == "-"):
            fh = sys.stdout
        else:
            fh = open(output,"w")
        
        i = 0
        
        # 1: write header
        fh.write("@HD	VN:1.0	SO:unsorted\n")
        for region in self.parse_regions():
            fh.write("@SQ	SN:"+region[0]+"	LN:"+str(region[2] - region[1] + 1)+"\n")
            
        fh.write("@PG	ID:0	PN:FlaiMapper_SSLM_to_SAM_conversion_script	VN:0.0\n")
        
        # 2: write alignment
        for region in self.get_alignment_files():
            logging.debug("   - Masked region: "+region[0])
            
            for read in self.parse_reads(region[1]):
                if(read.name):
                    fh.write(read.name)
                else:
                    fh.write("unknown_read_"+str(i))
                    i += 1
                
                strand = "60"
                fh.write("\t0\t"+region[0]+"\t"+str(read.start+1)+"\t"+strand+"\t"+str(read.stop - read.start)+"M\t*\t0\t0\t"+read.sequence+"\t*\tNH:i:1\n")
        
        fh.close()
