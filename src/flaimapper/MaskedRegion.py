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

import os,re,operator,argparse,sys


from flaimapper.Read import Read
from flaimapper.ncRNAfragment import ncRNAfragment


class MaskedRegion:
    """A masked region is a region masked in the reference genome to 
    indicate where ncRNAs are located.
    """
    def __init__(self,name,start,stop):
        self.name = name
        
        self.start = start
        self.stop = stop
    
    def reset(self):
        self.sequence = False
        
        self.start_positions = []
        self.stop_positions = []
    
    def parse_stats(self):
        self.reset()
        
        start_avg_lengths = []
        stop_avg_lengths = []
        
        for read in self.parse_reads():#Relies on inherented class, e.g. BAMParser or SSLMParser
            while(len(self.start_positions) < read.stop+1):				# Fix since 1.1.0: automatically scale  vector up if alignment falls outside range reference annotation
                self.start_positions.append(0)
                self.stop_positions.append(0)
                
                start_avg_lengths.append({})# Do an aggregated vector {21:10243} for 10243 observations of length 21
                stop_avg_lengths.append({})
            
            self.start_positions[read.start] += 1
            self.stop_positions[read.stop] += 1
            
            len_start = read.stop-read.start
            len_stop = read.start-read.stop
            
            if not start_avg_lengths[read.start].has_key(len_start):
                start_avg_lengths[read.start][len_start] = 0
            
            if not stop_avg_lengths[read.stop].has_key(len_stop):
                stop_avg_lengths[read.stop][len_stop] = 0
            
            start_avg_lengths[read.start][len_start] += 1
            stop_avg_lengths[read.stop][len_stop] += 1
        
        self.start_avg_lengths = []
        self.stop_avg_lengths = []
        
        for i in range(len(stop_avg_lengths)):
            avgLenF = self.get_median_of_map(start_avg_lengths[i])
            avgLenR = self.get_median_of_map(stop_avg_lengths[i])
            if(avgLenF):
                avgLenF = round(avgLenF+1)
            if(avgLenR):
                avgLenR = round(avgLenR-0.5)							# Why -0.5 -> because of rounding a negative number
            self.start_avg_lengths.append(avgLenF)
            self.stop_avg_lengths.append(avgLenR)
        
        return [self.start_positions,self.stop_positions,self.start_avg_lengths,self.stop_avg_lengths]
    
    def get_median_of_map(self,value_map):
        """
        input: 
        {
            21: 2
            23: 6
            24: 3
        }
        
        output:
            median of [21,21,23,23,23,23,23,23,24,24,24]
        
        ----
        
        trick:
        iter1:
        
        {
            21: 2  -- look at count of lowest key -> 2
            23: 6
            24: 3  -- look at count of larges key -> 3
        }
        
        remove the key with the lowest count and substract that number from the highest count:
        {   23: 6
            24: 1  }
        
        iter2:

        {   
            23: 6  -- key: 6
            24: 1  -- key: 1
         }
        
        remove lowest key:
        {   
            23: 5
         }
         
         median = 5
        
        
        ******
        
        in case of equal keys and n>2, remove both:
      { 21:5
        23:4
        24:5 }
        --> { 23:4 }
        
        ******
        in case of equal keys and n==2, return get_median([k1,k2])
        
      { 22:8
        26:8 }
        -> median([22,26]) == return (float(k1 + k2)) / 2
        """
        
        keys = sorted(value_map.keys())
        while len(keys) > 1:
            key_l = keys[0]
            key_t = keys[-1]
            
            if value_map[key_l] < value_map[key_t]:
                value_map[key_t] -= value_map[key_l]
                del(value_map[key_l])
                keys.remove(key_l)
            elif value_map[key_t] < value_map[key_l]:
                value_map[key_l] -= value_map[key_t]
                del(value_map[key_t])
                keys.remove(key_t)
            elif value_map[key_t] == value_map[key_l]:
                if len(keys) == 2:
                    return (float(keys[0] + keys[1])) / 2
                else:
                    del(value_map[key_l],value_map[key_t])
                    keys.remove(key_l)
                    keys.remove(key_t)
        
        if len(keys) > 0:
            return keys[0]
        else:
            return None
    
    def count_reads_per_region(self,fragments):							# @TODO change to 'sequencing_depth()'
        for fragment in fragments:
            fragment.supporting_reads = 0
        
        for read in self.parse_reads():
            for fragment in fragments:
                if(fragment.spans_read(read)):
                    fragment.add_supporting_reads(1)
