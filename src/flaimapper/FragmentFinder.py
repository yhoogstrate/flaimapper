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


import operator

from .ncRNAfragment import ncRNAfragment


class FragmentFinder:
    """Class Flaimapper Finds fragments in alignment files.
    """
    
    def __init__(self,masked_region,filter_parameters):
        #@todo Don't save this info indefinitely!!
        self.masked_region = masked_region
        self.filter_parameters = filter_parameters
        
        #@todo yield only, disable autorun, incorporate into bamparser or something to get:
        # for region in ...:
        #     for fragment in region.find_fragments():
        #         if output.format == gtf:
        #         fh.write(fragment.to_gtf())
        self.results = []
    
    def __iter__(self):
        for result in self.results:
            yield result
    
    def run(self):
        # Finds peaks
        peaksStart = self.find_peaks(self.masked_region.start_positions+[0])
        peaksStop = self.find_peaks(self.masked_region.stop_positions+[0])
        
        # Correct / filter noisy peaks
        peaksStart = self.smooth_filter_peaks(peaksStart)
        peaksStop = self.smooth_filter_peaks(peaksStop)
        
        # Trace start and stop positions together and obtain actual peaks
        self.results = self.find_fragments(
            peaksStart,
            peaksStop,
        
            self.masked_region.start_avg_lengths,
            self.masked_region.stop_avg_lengths)
    
    def find_peaks(self,plist,drop_cutoff=0.1):
        # Define variables:
        peaks = {}
        
        previous = 0
        highest = 0
        highestPos = -1
        
        # Walk over list of [start/stop]-position counts:
        for pos in range(len(plist)):
            current = plist[pos]
            if current > previous:# and (current > (noise_type_alpha_cutoff/100.0*max(plist)))):  
                if current > highest:
                    highest = current
                    highestPos = pos
            elif current < previous:
                #if (current < (drop_cutoff*highest)) and (highestPos != -1):
                #if (current < (100.0*drop_cutoff*highest)) and (highestPos != -1):   #
                #if (current < (10.0*highest)) and (highestPos != -1):   # 
                if (drop_cutoff * current < highest) and (highestPos != -1):
                    peaks[highestPos] = highest
                    
                    #highestPos = -1
                    highest = 0
            
            previous = current
        
        return peaks
    
    def smooth_filter_peaks(self,plist):
        """Smooth filtering
        """
        
        psorted = sorted(plist.iteritems(),key=operator.itemgetter(1),reverse=True)
        
        # There is a small mistake in the algorithm,
        # it should search not for ALL peaks
        # but only for ALL peaks except itself; position i can not be a noise product of i itself
        
        for i in range(len(psorted)):
            if(psorted[i] != False):
                item = psorted[i]
                for j in range(len(psorted)):							# Can be limited to size and -size of self.filter_parameters.matrix
                    if((psorted[j] != False) and (j != i)):
                        item2 = psorted[j]
                        diff = item2[0]-item[0]
                        if(self.filter_parameters.matrix.has_key(diff)):
                            perc = self.filter_parameters.matrix[diff]/100.0
                            if((perc*item[1]) > item2[1]):
                                psorted[j] = False
        
        pnew = {}
        
        for item in psorted:
            if(item != False):
                pnew[item[0]] = item[1]
        
        return pnew
    
    def find_fragments(self,pstart,pstop,pexpectedStart,pexpectedStop,prime_5_ext = 3,prime_3_ext=5,genomic_offset_masked_region=0):
        """Traceback:
        
        genomic_offset_masked_region - imagine your pre-miRNA is starts at position 400.000 in the genome; then your position should be 400.000 + start
        
        ----
        @return:
        @rtype:
        """
        fragments = []
        if(len(pstart) >= len(pstop)):									# More start than stop positions
            pstopSorted = sorted(pstop.iteritems(),key=operator.itemgetter(1))[::-1]
            for itema in pstopSorted:
                pos = itema[0]
                diff = pexpectedStop[pos]
                predictedPos = pos+diff+1								# 149 - 50 = 99; 149- 50 + 1 = 100 (example of read aligned to 100,149 (size=50)
                fragment = False
                
                highest = 0
                items = [s for s in pstart if ((s >= predictedPos-15) and (s <= predictedPos+15))]
                for item in items:
                    distance = abs(predictedPos - item)
                    penalty = 1.0 - (distance * 0.09)
                    score = pstart[item]*penalty 
                    if(score >= highest):
                        highest = pstart[item]
                        fragment = ncRNAfragment(item,pos,self.masked_region.region)
                        fragment.supporting_reads_start = pstart[fragment.start]
                        fragment.supporting_reads_stop = pstop[fragment.stop]
                
                if(fragment != False):
                    fragments.append(fragment)
                    del(pstart[fragment.start])
                    items = []
        else:															# More stop than start positions
            pstartSorted = sorted(pstart.iteritems(),key=operator.itemgetter(1))[::-1]
            for itema in pstartSorted:
                pos = itema[0]
                diff = pexpectedStart[pos]
                
                #@todo figure out if this requires << + 1
                predictedPos = pos+diff
                fragment = False
                
                highest = 0
                items = [s for s in pstop if ((s >= predictedPos-15) and (s <= predictedPos+15))]
                
                for item in items:
                    distance = abs(predictedPos - item)
                    penalty = 1.0 - (distance * 0.09)
                    score = pstop[item]*penalty 
                    if(score >= highest):
                        highest = pstop[item]
                        
                        fragment = ncRNAfragment(pos,item,self.masked_region.region)
                        fragment.supporting_reads_start = pstart[fragment.start]
                        fragment.supporting_reads_stop = pstop[fragment.stop]
                
                if(fragment != False):
                    fragments.append(fragment)
                    del(pstop[fragment.stop])
                    items = []
        #(counter >= fragment.start) and (counter < fragment.stop)
        
        return fragments
