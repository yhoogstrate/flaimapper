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



from FragmentContainer import FragmentContainer
from FragmentFinder import FragmentFinder
from AlignmentParser import AlignmentParser
from AlignmentDirectory import AlignmentDirectory



class FlaiMapperObject(FragmentContainer):
	def __init__(self,verbosity):
		self.verbosity = verbosity
		
		self.alignment_directories = []
		self.alignment_directories_indexed = {}
		
		self.sequences = {}
		
		if(self.verbosity == "verbose"):
			print " - Initiated FlaiMapper Object"
	
	def add_alignment_directory(self,alignment_directory):
		if(alignment_directory.__class__.__name__ != AlignmentDirectory.__name__):
			raise TypeError, "alignment_directory is not of class type " + AlignmentDirectory.__name__ + " but of class type: " + alignment_directory.__class__.__name__
		else:
			self.alignment_directories.append(alignment_directory)
	
	def index(self):
		if(self.verbosity == "verbose"):
			print " - Indexing alignment files"
		
		for alignment_directory in self.alignment_directories:
			if(self.verbosity == "verbose"):
				print "   - Indexing "+alignment_directory.path
			for ncRNA in alignment_directory.index.keys():
				if(not self.alignment_directories_indexed.has_key(ncRNA)):
					self.alignment_directories_indexed[ncRNA] = []
				self.alignment_directories_indexed[ncRNA].append(alignment_directory.index[ncRNA])
	
	def run(self):
		self.index()
		
		if(self.verbosity == "verbose"):
			print " - Running fragment detection"
		
		for ncRNA in self.alignment_directories_indexed.keys():
			if(self.verbosity == "verbose"):
				print "   - Analysing: "+ncRNA
			
			aligned_reads = AlignmentParser(ncRNA,self.alignment_directories_indexed[ncRNA],self.verbosity)
			aligned_reads.parse_stats()
			
			predicted_fragments = FragmentFinder(ncRNA,aligned_reads)#,True)
			predicted_fragments.run()
			
			self.add_fragments(predicted_fragments)
	
	def count_reads_per_region_custom_table(self,regions,links,all_predicted_fragments,reference_offset=0):
		"""
		All sequences in our library of ncRNAs have been extended with 10 bases.
		"""
		
		stats_table = {}
		stats_table['experimental']     = {'error_5p':{"<-5":0,-5:0,-4:0,-3:0,-2:0,-1:0,0:0,1:0,2:0,3:0,4:0,5:0,">5":0},'error_3p':{"<-5":0,-5:0,-4:0,-3:0,-2:0,-1:0,0:0,1:0,2:0,3:0,4:0,5:0,">5":0},'predicted':0,'not_predicted_no_reads':0,'not_predicted_with_reads':0}
		stats_table['not_experimental'] = {'error_5p':{"<-5":0,-5:0,-4:0,-3:0,-2:0,-1:0,0:0,1:0,2:0,3:0,4:0,5:0,">5":0},'error_3p':{"<-5":0,-5:0,-4:0,-3:0,-2:0,-1:0,0:0,1:0,2:0,3:0,4:0,5:0,">5":0},'predicted':0,'not_predicted_no_reads':0,'not_predicted_with_reads':0}
		
		
		#a = c(1:10) mse_a = sum((a - mean(a)) ^ 2) / length(a)
		
		if(self.verbosity == "verbose"):
			print " - Running fragment detection"
		
		i = 0
		j = 0
		
		
		for ncRNA in all_predicted_fragments.keys():
			if(links.has_key(ncRNA)):
				if(self.verbosity == "verbose"):
					print "   - Analysing: "+ncRNA
				
				annotations = regions.index[links[ncRNA]]
				
				predicted_fragments = all_predicted_fragments[ncRNA].getResults()
				
				#aligned_reads = AlignmentParser(ncRNA,self.alignment_directories_indexed[ncRNA],self.verbosity)
				#aligned_reads.count_reads_per_region(annotations)
				#aligned_reads.parse_stats()
				
				#predicted_fragments_obj = FragmentFinder(ncRNA,aligned_reads)
				#predicted_fragments_obj.run()
				#predicted_fragments = predicted_fragments_obj.getResults()
				
				i += 1
				
				for annotation in annotations.fragments:
					closest = self.find_closest_overlapping_fragment(annotation,predicted_fragments,reference_offset)
					j += 1
					
					if(closest):
						errors = self.find_errors(annotation,closest)
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
		
		print i,"annotated pre-miRNAs"
		print j,"annotated miRNAs"
		
		return stats_table
	
	def count_reads_per_region_custom_mse(self,regions,links,all_predicted_fragments,reference_offset=0):
		"""
		All sequences in our library of ncRNAs have been extended with 10 bases.
		"""
		
		
		err_5p = []
		err_3p = []
		
		
		if(self.verbosity == "verbose"):
			print " - Running fragment detection"
		
		i = 0
		j = 0
		
		import numpy
		
		for ncRNA in all_predicted_fragments.keys():
			if(links.has_key(ncRNA)):
				if(self.verbosity == "verbose"):
					print "   - Analysing: "+ncRNA
				
				match = re.search("chr[^:]+:([0-9]+)-([0-9]+):",ncRNA)
				seq_length = abs(int(match.group(1)) - int(match.group(2)))
				
				annotations = regions.index[links[ncRNA]]
				
				predicted_fragments = all_predicted_fragments[ncRNA].getResults()
				
				
				i += 1
				
				for annotation in annotations.fragments:
					closest = self.find_closest_overlapping_fragment(annotation,predicted_fragments,reference_offset)
					j += 1
					
					if(closest):
						errors = self.find_errors(annotation,closest)
						err_5p.append(errors[0])
						err_3p.append(errors[1])
						
					else:
						ann_length = abs(annotation.stop - annotation.start) * 2
						#err_5p.append(ann_length)# take missed miRNA length instead
						#err_3p.append(ann_length)# take missed miRNA length instead
						err_5p.append(seq_length)# take missed miRNA length instead
						err_3p.append(seq_length)# take missed miRNA length instead
		
		# Calc Root Mean Square Error (RMSQ)
		err_5p = numpy.sqrt(numpy.mean(numpy.array(err_5p)**2))
		err_3p = numpy.sqrt(numpy.mean(numpy.array(err_3p)**2))
		
		return [err_5p,err_3p]
	
	def count_reads_per_region(self,regions,links,reference_offset=0):
		"""
		All sequences in our library of ncRNAs have been extended with 10 bases.
		"""
		self.index()
		
		stats_table = {}
		stats_table['experimental']     = {'error_5p':{"<-5":0,-5:0,-4:0,-3:0,-2:0,-1:0,0:0,1:0,2:0,3:0,4:0,5:0,">5":0},'error_3p':{"<-5":0,-5:0,-4:0,-3:0,-2:0,-1:0,0:0,1:0,2:0,3:0,4:0,5:0,">5":0},'predicted':0,'not_predicted_no_reads':0,'not_predicted_with_reads':0}
		stats_table['not_experimental'] = {'error_5p':{"<-5":0,-5:0,-4:0,-3:0,-2:0,-1:0,0:0,1:0,2:0,3:0,4:0,5:0,">5":0},'error_3p':{"<-5":0,-5:0,-4:0,-3:0,-2:0,-1:0,0:0,1:0,2:0,3:0,4:0,5:0,">5":0},'predicted':0,'not_predicted_no_reads':0,'not_predicted_with_reads':0}
		
		if(self.verbosity == "verbose"):
			print " - Running fragment detection"
		
		i = 0
		j = 0
		
		for ncRNA in self.alignment_directories_indexed.keys():
			if(links.has_key(ncRNA)):
				if(self.verbosity == "verbose"):
					print "   - Analysing: "+ncRNA
				
				annotations = regions.index[links[ncRNA]]
				
				aligned_reads = AlignmentParser(ncRNA,self.alignment_directories_indexed[ncRNA],self.verbosity)
				aligned_reads.count_reads_per_region(annotations)
				aligned_reads.parse_stats()
				
				predicted_fragments_obj = FragmentFinder(ncRNA,aligned_reads)
				predicted_fragments_obj.run()
				predicted_fragments = predicted_fragments_obj.getResults()
				
				i += 1
				
				for annotation in annotations.fragments:
					closest = self.find_closest_overlapping_fragment(annotation,predicted_fragments,reference_offset)
					j += 1
					
					if(closest):
						errors = self.find_errors(annotation,closest)
						err_5p = errors[0]
						err_3p = errors[1]
						
						#if(abs(err_5p) > 5):
						#	print err_5p
						
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
		
		print i,"annotated pre-miRNAs"
		print j,"annotated miRNAs"
		
		return stats_table
	
	def count_error_with_intensity(self,regions,links,reference_offset=0):
		"""
		All sequences in our library of ncRNAs have been extended with 10 bases.
		"""
		self.index()
		
		out = []
		
		if(self.verbosity == "verbose"):
			print " - Running fragment detection"
		
		for ncRNA in self.alignment_directories_indexed.keys():
			if(links.has_key(ncRNA)):
				if(self.verbosity == "verbose"):
					print "   - Analysing: "+ncRNA
				
				annotations = regions.index[links[ncRNA]]
				
				aligned_reads = AlignmentParser(ncRNA,self.alignment_directories_indexed[ncRNA],self.verbosity)
				aligned_reads.count_reads_per_region(annotations)
				aligned_reads.parse_stats()
				
				predicted_fragments_obj = FragmentFinder(ncRNA,aligned_reads)
				predicted_fragments_obj.run()
				predicted_fragments = predicted_fragments_obj.getResults()
				
				for annotation in annotations.fragments:
					closest = self.find_closest_overlapping_fragment(annotation,predicted_fragments,reference_offset)
					
					if(closest):
						errors = self.find_errors(annotation,closest)
						err_5p = errors[0]
						err_3p = errors[1]
						
						out.append({'5p':[closest[2],err_5p],'3p':[closest[3],err_3p]})
		
		return out
	
	def find_closest_overlapping_fragment(self,annotated_fragment,predicted_fragments,reference_offset=0):
		closest = False
		closest_overlapping_bases = 0
		
		for predicted_fragment in predicted_fragments:
			annotated = [annotated_fragment.start,annotated_fragment.stop]
			if(predicted_fragment.has_key('start_supporting_reads')):
				predicted = [(predicted_fragment["start"] - reference_offset), (predicted_fragment["stop"] - reference_offset),predicted_fragment['start_supporting_reads'],predicted_fragment['stop_supporting_reads']]
			else:
				predicted = [(predicted_fragment["start"] - reference_offset), (predicted_fragment["stop"] - reference_offset),0,0]
			
			overlap = self.find_overlapping_bases(annotated,predicted)
			if(overlap > 0 and overlap > closest_overlapping_bases):
				closest_overlapping_bases = overlap
				closest = predicted
		
		return closest
	
	def find_overlapping_bases(self,fragment_1,fragment_2):
		if(fragment_2[0] < fragment_1[0]):
			return self.find_overlapping_bases(fragment_2,fragment_1)
		else:
			return fragment_1[1] - fragment_2[0]
	
	"""
	Example:
		   [ miRNA ]
		[ fragment* ]
	
	mirna: 4,12
	fragment: 0,14
	
	error_5p = 0 - 4 = -4
	error_3p = 13 - 12 = 1
	
	"""
	def find_errors(self,annotated_fragment,predicted_fragment):
		
		error_5p = predicted_fragment[0] - annotated_fragment.start
		error_3p = predicted_fragment[1] - annotated_fragment.stop
		
		return [error_5p,error_3p]
