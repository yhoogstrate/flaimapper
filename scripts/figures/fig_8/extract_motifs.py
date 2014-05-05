#!/bin/env python

"""
This script finds the suffixes of the fragments predicted on H/ACA-box
snoRNAs located on the 3' halve of the H/ACA-box snoRNA.
"""

# Offset how far a fragment must relatively be compared to the ncRNAs length
relative_position_boundary = 0.6


# Read the ncRNA references

def get_lengths():
	fh = open("../../../share/annotations/ncRNA_annotation/ncRNdb09_with_tRNAs_and_Pseudogenes__21_oct_2011__hg19.fasta","r")
	lines = fh.read().split("\n")
	fh.close()
	
	sequence_lengths = {}
	
	for line in lines:
		line = line.strip()
		if(len(line) > 0):
			if(line[0] == ">"):
				name = line
				sequence_lengths[name] = 0
			else:
				sequence_lengths[name] += len(line)
	
	return sequence_lengths

sequence_lengths = get_lengths()


def complete_sequence(sequence,length,direction,fill="N"):
	diff = max(0,length - len(sequence))
	fill = diff * fill
	sequence = sequence.upper().replace('T','U')
	if(direction == "suffix"):
		return fill + str(sequence)
	elif(direction == "prefix"):
		return str(sequence) + fill



def get_sequences(arg_subtype):
	# Find corresponding fragments and find suffixes
	fh = open("../../../output/FlaiMapper/SRP002175/02_output_flaimapper.txt","r")
	lines = fh.read().split("\n")
	fh.close()
	
	sequences = []
	max_length = 0
	
	for line in lines:
		line = line.strip()
		if(len(line) > 0):
			params = line.split("\t")
			
			if params[1].find(arg_subtype) > -1:
				s = {'name':params[1],'sequence':params[4],'length': len(params[4]),'precursor_length':sequence_lengths[params[1]],'abs_center':(float(params[2])+float(params[3])) * 0.5}
				s['relative_position'] = s['abs_center'] / s['precursor_length']
				sequences.append(s)
				if(s['length'] > max_length):
					max_length = s['length']
	
	return max_length, sequences



for subtype in ["MIR", "SNORA", "SNORD"]:
	print "Analysing subtype: "+subtype
	max_length, sequences = get_sequences(subtype)
	for direction in ["prefix", "suffix"]:
		fh = open("motifs_"+subtype+"_"+direction+".txt","w")
		print " - Calculating "+direction+"es"
		for sequence in sequences:
			completed_sequence = complete_sequence(sequence['sequence'],max_length,direction,"N")
			if(direction == "prefix" and sequence['relative_position'] <= relative_position_boundary):
				fh.write(completed_sequence[:18]+"\n")
			elif(direction == "suffix" and sequence['relative_position'] >= relative_position_boundary):
				fh.write(completed_sequence[-18:]+"\n")
		fh.close()

