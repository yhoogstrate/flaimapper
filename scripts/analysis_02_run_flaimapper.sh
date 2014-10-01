#!/bin/bash

echo "02. Annotating fragments"


mkdir -p ../output/FlaiMapper/SRP002175_merged
mkdir -p ../output/FlaiMapper/SRP002175
mkdir -p ../output/FlaiMapper/SRP006788
mkdir -p ../output/FlaiMapper/SRP028959									# 4-1x Ion Torrent PGM HeLa merged
mkdir -p ../output/FlaiMapper/SRP034013_merged							# 3x HiSeq2000 B-cells merged
mkdir -p ../output/FlaiMapper/SRP041082_merged							# 2x HiSeq2000 PCa merged

# ---------------------------SRP002175----------------------------------

echo "  - Dataset: SRP002175 (experiments merged)"

flaimapper \
	 -v \
	 -o ../output/FlaiMapper/SRP002175_merged/01_output_flaimapper.txt \
	 -m ../share/annotations/ncRNA_annotation/ncrnadb09.gtf \
	 -r ../share/annotations/ncRNA_annotation/ncrnadb09.fa \
	 -f 1 \
		../share/small_RNA-seq_alignments/SRP002175/SRR038852.bam \
		../share/small_RNA-seq_alignments/SRP002175/SRR038853.bam \
		../share/small_RNA-seq_alignments/SRP002175/SRR038854.bam \
		../share/small_RNA-seq_alignments/SRP002175/SRR038855.bam \
		../share/small_RNA-seq_alignments/SRP002175/SRR038856.bam \
		../share/small_RNA-seq_alignments/SRP002175/SRR038857.bam \
		../share/small_RNA-seq_alignments/SRP002175/SRR038858.bam \
		../share/small_RNA-seq_alignments/SRP002175/SRR038859.bam \
		../share/small_RNA-seq_alignments/SRP002175/SRR038860.bam \
		../share/small_RNA-seq_alignments/SRP002175/SRR038861.bam \
		../share/small_RNA-seq_alignments/SRP002175/SRR038862.bam \
		../share/small_RNA-seq_alignments/SRP002175/SRR038863.bam

# ----------------------------------------------------------------------

echo "  - Dataset: SRP002175"

SAMPLES=( "SRR038852" "SRR038853" "SRR038854" "SRR038855" "SRR038856" "SRR038857" "SRR038858" "SRR038859" "SRR038860" "SRR038861" "SRR038862" "SRR038863" )

for SAMPLE in "${SAMPLES[@]}"
do :
	flaimapper \
		 -v \
		 -o "../output/FlaiMapper/SRP002175/01_output_flaimapper_"$SAMPLE".txt" \
		 -m ../share/annotations/ncRNA_annotation/ncrnadb09.gtf \
		 -r ../share/annotations/ncRNA_annotation/ncrnadb09.fa \
		 -f 1 \
			"../share/small_RNA-seq_alignments/SRP002175/"$SAMPLE".bam"
done

unset SAMPLE SAMPLES

# ----------------------------------------------------------------------

# ---------------------------SRP006788----------------------------------

echo "  - Dataset: SRP006788"

SAMPLES=( "SRR207111_HeLa18-30" "SRR207112_HeLa18-30_RRP40" "SRR207113_HeLa18-30_AGO1_2" "SRR207114_HeLa18-30_AGO1_2_RRP40" "SRR207115_HeLa18-30_XRN1_2" "SRR207116_HeLa18-30_N" )

for SAMPLE in "${SAMPLES[@]}"
do :
	flaimapper \
		 -v \
		 -o "../output/FlaiMapper/SRP006788/01_output_flaimapper_"$SAMPLE".txt" \
		 -m ../share/annotations/ncRNA_annotation/ncrnadb09.gtf \
		 -r ../share/annotations/ncRNA_annotation/ncrnadb09.fa \
		 -f 1 \
			"../share/small_RNA-seq_alignments/SRP006788/"$SAMPLE".bam"
done

unset SAMPLE SAMPLES

# ----------------------------------------------------------------------
