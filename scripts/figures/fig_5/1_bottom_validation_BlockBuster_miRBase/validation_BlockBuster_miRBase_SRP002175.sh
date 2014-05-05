#!/bin/bash

# Step 01: annotate the fragments
../../../../src/sslm2bed-converter.py \
                  -o ../../../../share/small_RNA-seq_alignments/SRP002175/blockbuster_input_SRP002175_merged.unsorted.bed \
                  ../../../../share/small_RNA-seq_alignments/SRP002175/SRR038853 \
                  ../../../../share/small_RNA-seq_alignments/SRP002175/SRR038855 \
                  ../../../../share/small_RNA-seq_alignments/SRP002175/SRR038857 \
                  ../../../../share/small_RNA-seq_alignments/SRP002175/SRR038859 \
                  ../../../../share/small_RNA-seq_alignments/SRP002175/SRR038861 \
                  ../../../../share/small_RNA-seq_alignments/SRP002175/SRR038863 \
                  ../../../../share/small_RNA-seq_alignments/SRP002175/SRR038852 \
                  ../../../../share/small_RNA-seq_alignments/SRP002175/SRR038854 \
                  ../../../../share/small_RNA-seq_alignments/SRP002175/SRR038856 \
                  ../../../../share/small_RNA-seq_alignments/SRP002175/SRR038858 \
                  ../../../../share/small_RNA-seq_alignments/SRP002175/SRR038860 \
                  ../../../../share/small_RNA-seq_alignments/SRP002175/SRR038862

sort -k 1,1 -k 2,2n -k 3,3n ../../../../share/small_RNA-seq_alignments/SRP002175/blockbuster_input_SRP002175_merged.unsorted.bed > \
                            ../../../../share/small_RNA-seq_alignments/SRP002175/blockbuster_input_SRP002175_merged.sorted.bed
rm ../../../../share/small_RNA-seq_alignments/SRP002175/blockbuster_input_SRP002175_merged.unsorted.bed


# Step 02: alignment for several parameters:
for sd in {0.05,0.1,0.2,0.35,0.5,1,1.5,2}
do
	for dist in {1..50}
	do
		echo "blockbuster_output.min-dist_"$dist"_sd_"$sd".txt"
		
		blockbuster.x \
			-format 1 \
			-scale $sd \
			-distance $dist \
			"../../../../share/small_RNA-seq_alignments/SRP002175/blockbuster_input_SRP002175_merged.sorted.bed" \
			> "../../../../output/BlockBuster/SRP002175/min-dist_"$dist"_sd_"$sd".txt"
	done
done
