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
				
				start_avg_lengths.append([])
				stop_avg_lengths.append([])
			
			self.start_positions[read.start] += 1
			self.stop_positions[read.stop] += 1
			
			start_avg_lengths[read.start].append(read.stop-read.start)
			stop_avg_lengths[read.stop].append(read.start-read.stop)
		
		self.start_avg_lengths = []
		self.stop_avg_lengths = []
		
		for i in range(len(stop_avg_lengths)):
			avgLenF = self.get_median(start_avg_lengths[i])
			avgLenR = self.get_median(stop_avg_lengths[i])
			if(avgLenF):
				avgLenF = round(avgLenF+1)
			if(avgLenR):
				avgLenR = round(avgLenR-0.5)							# Why -0.5 -> because of rounding a negative number
			self.start_avg_lengths.append(avgLenF)
			self.stop_avg_lengths.append(avgLenR)
		
		return [self.start_positions,self.stop_positions,self.start_avg_lengths,self.stop_avg_lengths]
	
	def parse_reads_stacked(self,return_sorted = True):
		if(return_sorted):
			index = {}
			for read in self.parse_reads():
				if(not read.start in index.keys()):
					index[read.start] = {}
				
				if(read.stop in index[read.start].keys()):
					index[read.start][read.stop] += 1
				else:
					index[read.start][read.stop] = 1
			
			for start in sorted(index.keys()):
				for stop in sorted(index[start].keys()):
					yield [Read(start,stop,None),index[start][stop]]
		else:
			index = {}
			for read in self.parse_reads():
				if(not read.start in index.keys()):
					index[read.start] = {}
				
				if(read.stop in index[read.start].keys()):
					index[read.start][read.stop] += 1
				else:
					index[read.start][read.stop] = 1
			
			for start in index.keys():
				for stop in index[start].keys():
					yield [Read(start,stop,None),index[start][stop]]
	
	def get_median(self,numericValues):
		"""Finds the median of a vector"""								# @TODO move to utils and rename to 'median()'
		theValues = sorted(numericValues)
		count = len(theValues)
		
		if count == 0:
			return None
		elif count % 2 == 1:
			return theValues[(count+1)/2-1]
		else:
			lower = theValues[count/2-1]
			upper = theValues[count/2]
			return (float(lower + upper)) / 2
	
	def count_reads_per_region(self,fragments):							# @TODO change to 'sequencing_depth()'
		for fragment in fragments:
			fragment.supporting_reads = 0
		
		for read in self.parse_reads():
			for fragment in fragments:
				if(fragment.spans_read(read)):
					fragment.add_supporting_reads(1)
