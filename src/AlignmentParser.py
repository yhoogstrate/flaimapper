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
 - MSc. Youri Hoogstrate
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


__version_info__ = ('1', '0', '0')
__version__ = '.'.join(__version_info__)
__author__ = 'Youri Hoogstrate'
__homepage__ = 'https://github.com/yhoogstrate/flaimapper'
__license__ = 'GPL3'



import os,re,random,operator,argparse,sys



class AlignmentParser:
	"""parseNcRNA is a class that parses the SSLM alignment files.
	"""
	parseNcRNA_prog = re.compile('.*?LOCI=\[(chr[^:]+):([0-9]+)-([0-9]+):')# Shouldn't this be in the class self.prog?
		
	def __init__(self,ncRNA,fileList,verbosity):
		self.verbosity = verbosity
		self.positions = re.match(self.parseNcRNA_prog,ncRNA.replace('location=',''))
		
		self.name = ncRNA
		self.file_list = fileList
	
	def reset(self):
		self.sequence = False
		
		self.start_avg_lengths = []
		self.stop_avg_lengths = []
		
		self.coverage = []
	
	def parse_coverage(self):
		self.reset()
		"""parse
		
		----
		@param self.name: The name of the GENE, neccesairy for finding the gene-size easily.
		@param self.file_list: 
		"""
		
		sequence_length = int(self.positions.group(3))-int(self.positions.group(2))+1
		
		fileHandles = []
		fileContents = []
		
		for file in self.file_list:
			fh = open(file,'r')
			fileHandles.append(fh)
			contents = fh.readlines()
			if(not self.sequence):
				self.sequence = contents[1].strip()
				if(sequence_length != len(self.sequence) and self.verbosity != "quiet" and self.name.upper().find("TRNA") == -1):# tRNA annotations have alternative splicing and do thus not match
					print "***Warning: (ncRNA: "+self.name+") sequence length doesn't match sequence: "+str(sequence_length)+" != "+str(len(self.sequence))
					return False
			fileContents.append(contents[2:])
		
		for i in range(len(self.sequence)):
			self.coverage.append(0)
		
		for file in fileContents:
			for i in range(len(file)/2):
				numberofhits = int(file[i*2].strip()[::-1].split('STIH',2)[0][::-1].replace('revcomp',''))
				read = file[(i*2)+1].strip()
				
				startPos = self.getStartPos(read)
				stopPos = self.getStopPos(read)
				
				if(stopPos > len(self.sequence)):
					while(len(self.coverage) < stopPos):
						self.coverage.append(0)
					
					self.sequence_length = stopPos
				
				for i in range(startPos,stopPos):
					self.coverage[i] += 1
		
		for fh in fileHandles:
			fh.close()
	
	def parse_reads(self):
		self.reset()
		"""parse
		
		----
		@param self.name: The name of the GENE, neccesairy for finding the gene-size easily.
		@param self.file_list: 
		"""
		
		reads = {}
		
		sequence_length = int(self.positions.group(3))-int(self.positions.group(2))+1
		
		fileHandles = []
		fileContents = []
		
		for file in self.file_list:
			fh = open(file,'r')
			fileHandles.append(fh)
			contents = fh.readlines()
			if(not self.sequence):
				self.sequence = contents[1].strip()
				if(sequence_length != len(self.sequence) and self.verbosity != "quiet" and self.name.upper().find("TRNA") == -1):# tRNA annotations have alternative splicing and do thus not match
					print "***Warning: (ncRNA: "+self.name+") sequence length doesn't match sequence: "+str(sequence_length)+" != "+str(len(self.sequence))
					return False
			fileContents.append(contents[2:])
		
		for file in fileContents:
			for i in range(len(file)/2):
				numberofhits = int(file[i*2].strip()[::-1].split('STIH',2)[0][::-1].replace('revcomp',''))
				read = file[(i*2)+1].strip()
				
				startPos = self.getStartPos(read)
				stopPos = self.getStopPos(read)
				
				read_id = str(startPos)+"_"+str(stopPos)
				
				if(not reads.has_key(read_id)):
					reads[read_id] = {'startPos':startPos,'hits':0,'stopPos':stopPos}
				reads[read_id]['hits'] += numberofhits
		
		for fh in fileHandles:
			fh.close()
		
		return reads
	
	def parse_stats(self):
		self.reset()
		
		sequence_length = int(self.positions.group(3))-int(self.positions.group(2))+1
		
		self.start_positions = []
		self.stop_positions = []
		
		for i in range(sequence_length+1):
			self.start_positions.append(0)
			self.stop_positions.append(0)
			self.start_avg_lengths.append([])
			self.stop_avg_lengths.append([])
		
		fileHandles = []
		fileContents = []
		
		for file in self.file_list:
			fh = open(file,'r')
			fileHandles.append(fh)
			contents = fh.readlines()
			if(not self.sequence):
				self.sequence = contents[1].strip()
				if(sequence_length != len(self.sequence) and self.verbosity != "quiet" and self.name.upper().find("TRNA") == -1):# tRNA annotations have alternative splicing and do thus not match
					print "***Warning: (ncRNA: "+self.name+") sequence length doesn't match sequence: "+str(sequence_length)+" != "+str(len(self.sequence))
					
			fileContents.append(contents[2:])
		
		for file in fileContents:
			for i in range(len(file)/2):
				
				numberofhits = int(file[i*2].strip()[::-1].split('STIH',2)[0][::-1].replace('revcomp',''))
				read = file[(i*2)+1].strip()
				
				startPos = self.getStartPos(read)
				stopPos = self.getStopPos(read)
				if(stopPos > sequence_length):
					stopPos = sequence_length
				
				self.start_positions[startPos] += numberofhits
				self.stop_positions[stopPos] += numberofhits
				self.start_avg_lengths[startPos].append(stopPos-startPos)
				self.stop_avg_lengths[stopPos].append(startPos-stopPos)
				
		
		for i in range(len(self.start_avg_lengths)):
			avgLenF = self.getMedian(self.start_avg_lengths[i])
			avgLenR = self.getMedian(self.stop_avg_lengths[i])
			if(avgLenF):
				avgLenF = round(avgLenF)
			if(avgLenR):
				avgLenR = round(avgLenR-0.5)# Why -0.5 -> because of rounding
			self.start_avg_lengths[i] = avgLenF
			self.stop_avg_lengths[i] = avgLenR
		
		for fh in fileHandles:
			fh.close()
	
	def get_start_positions(self):
		return self.start_positions
		
	def get_stop_postions(self):
		return self.stop_positions
	
	def get_start_avg_lengths(self):
		return self.start_avg_lengths
	
	def get_stop_avg_lengths(self):
		return self.stop_avg_lengths
	
	
	def getMedian(self,numericValues):
		"""Finds the median of a vector"""
		theValues = sorted(numericValues)
		count = len(theValues)
		if count == 0:
			return False
		elif count % 2 == 1:
			return theValues[(count+1)/2-1]
		else:
			lower = theValues[count/2-1]
			upper = theValues[count/2]
			return (float(lower + upper)) / 2
	
	def getStartPos(self,read,extention='-'):
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
	
	def getStopPos(self,read,extention='-'):
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
	
	def get_sequence(self):
		return self.sequence
	
	def getResults(self):
		"""
		Format:
		
		{'startPositions':startPositions,'stopPositions':stopPositions,'startAvgLengths':lengthsF,'stopAvgLengths':lengthsR}
		
		----
		@return:
		@rtype: dictrionary
		"""
		return self.results
	
	def count_reads_per_region(self,annotations):
		self.reset()
		
		sequence_length = int(self.positions.group(3))-int(self.positions.group(2))+1
		
		self.start_positions = []
		self.stop_positions = []
		
		fileHandles = []
		fileContents = []
		
		for file in self.file_list:
			fh = open(file,'r')
			fileHandles.append(fh)
			contents = fh.readlines()
			if(not self.sequence):
				self.sequence = contents[1].strip()
				if(sequence_length != len(self.sequence) and self.verbosity != "quiet" and self.name.upper().find("TRNA") == -1):# tRNA annotations have alternative splicing and do thus not match
					print "***Warning: (ncRNA: "+self.name+") sequence length doesn't match sequence: "+str(sequence_length)+" != "+str(len(self.sequence))
					
			fileContents.append(contents[2:])
		
		for file in fileContents:
			for i in range(len(file)/2):
				numberofhits = int(file[i*2].strip()[::-1].split('STIH',2)[0][::-1].replace('revcomp',''))
				read = file[(i*2)+1].strip()
				
				startPos = self.getStartPos(read)
				stopPos = self.getStopPos(read)
				if(stopPos > sequence_length):
					stopPos = sequence_length
				
				for annotation in annotations.fragments:
					if(self.has_overlap([annotation.start,annotation.stop],[startPos,stopPos])):
						annotation.add_supporting_reads(numberofhits)
		
		for fh in fileHandles:
			fh.close()
	
	def has_overlap(self,left_obj,right_obj):
		"""
		Sort the objects on the first pos
		
		[     ]
		    [     ]
		
		[         ]
		    [    ]
		
		Then see if the start position of the second falls in the range of the first.
		"""
		
		if(right_obj[0] < left_obj[0]):
			return self.has_overlap(right_obj,left_obj)
		else:
			return (right_obj[0] >= left_obj[0] and right_obj[0] <= left_obj[1])
