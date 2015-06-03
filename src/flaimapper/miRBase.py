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

from flaimapper.ncRNA import ncRNA
from flaimapper.ncRNAfragment import ncRNAfragment

class miRBase:
	def __init__(self,arg_filename):
		self.mirs = []
		self.index = {}
		self.parse(arg_filename)
	
	def parse(self,arg_filename):
		fh = open(arg_filename,'r')
		content = fh.readlines()
		fh.close()
		
		state = 'closed'
		
		for line in content:
			line = line.rstrip()
			key = line[0:5].rstrip()
			if(state == 'closed'):
				if(key == 'DE'):
					organism = ' '.join(line[5:].split(' ',3)[0:2])
					if(organism == 'Homo sapiens'):
						info = {}
						info['name'] = line.split('Homo sapiens',1)[1].replace('stem-loop','').strip()
						info['mirs'] = []
						info['seq'] = ''
						info['aliases'] = []
						info['name_mir'] = 'unknown'
						state = 'open'
			elif(state == 'open'):
					if(key == "DR"):
						ids = []
						tmp_ids = line[5:].split(";")
						
						for i in range(len(tmp_ids)):
							tmp_ids[i] = tmp_ids[i].strip()
						
						if(tmp_ids[0] == "ENTREZGENE"):
							for tmp_id in tmp_ids[1:]:
								ids.append(tmp_id.rstrip("."))
						
						info['aliases'] = ids
					
					if(key == 'FT'):
						sline = line[5:].lstrip()
						if((sline[0:5] == 'miRNA') or (sline[0:8] == 'misc_RNA')):
							info['positions'] = sline.replace('miRNA','').replace('misc_RNA','').replace(' ','').split("..")
							info['positions'] = {'start':int(info['positions'][0])-1,'stop':int(info['positions'][1])}
						elif(sline[0:len('/product=')] == '/product='):
							info['name_mir'] = sline.split('=',1)[1].strip('"')
						elif(sline[0:len('/evidence=')] == '/evidence='):
							evidence = sline.split("=",1)[1]
							info['mirs'].append({'name':info['name'],'pos':info['positions'],'evidence':evidence})
							del(info['positions'])
							del(info['name_mir'])
					elif(key == ''):
						info['seq'] += line[5:].replace(' ','').replace('1','').replace('2','').replace('3','').replace('4','').replace('5','').replace('6','').replace('7','').replace('8','').replace('9','').replace('0','').upper().replace('U','T')
					elif(key == '//'):
						state = 'closed'
						ncRNAObj = ncRNA(info['name'])
						for fragment in info['mirs']:
							fragmentObj = ncRNAfragment(fragment['pos']['start'],fragment['pos']['stop'],None,None)
							fragmentObj.set_sequence(info['seq'][fragment['pos']['start']:fragment['pos']['stop']])
							fragmentObj.set_name(fragment['name'])
							fragmentObj.set_evidence(fragment['evidence'])
							ncRNAObj.add_fragments(fragmentObj)
						ncRNAObj.set_parameter("aliases",info["aliases"])
						self.mirs.append(ncRNAObj)
						self.index[info['name']] = ncRNAObj
						info = {}
	
	def get_miRNAs(self):
		return self.mirs
