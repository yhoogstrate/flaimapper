#!/bin/bash

# Step 01: annotate the fragments
../src/flaimapper.py \
                     -v \
                     -o ../output/SRP002175/01_output_flaimapper.txt \
                     -f 1 \
                        ../share/small_RNA-seq_alignments/SRP002175/SRR038853 \
                        ../share/small_RNA-seq_alignments/SRP002175/SRR038855 \
                        ../share/small_RNA-seq_alignments/SRP002175/SRR038857 \
                        ../share/small_RNA-seq_alignments/SRP002175/SRR038859 \
                        ../share/small_RNA-seq_alignments/SRP002175/SRR038861 \
                        ../share/small_RNA-seq_alignments/SRP002175/SRR038863 \
                        ../share/small_RNA-seq_alignments/SRP002175/SRR038852 \
                        ../share/small_RNA-seq_alignments/SRP002175/SRR038854 \
                        ../share/small_RNA-seq_alignments/SRP002175/SRR038856 \
                        ../share/small_RNA-seq_alignments/SRP002175/SRR038858 \
                        ../share/small_RNA-seq_alignments/SRP002175/SRR038860 \
                        ../share/small_RNA-seq_alignments/SRP002175/SRR038862

../src/flaimapper.py \
                     -v \
                     -o ../output/SRP006788/01.a_output_flaimapper.txt \
                     -f 1 \
                        ../share/small_RNA-seq_alignments/SRP006788/SRR207111_HeLa18-30

../src/flaimapper.py \
                     -v \
                     -o ../output/SRP006788/01.b_output_flaimapper.txt \
                     -f 1 \
                        ../share/small_RNA-seq_alignments/SRP006788/SRR207112_HeLa18-30_RRP40

../src/flaimapper.py \
                     -v \
                     -o ../output/SRP006788/01.c_output_flaimapper.txt \
                     -f 1 \
                        ../share/small_RNA-seq_alignments/SRP006788/SRR207113_HeLa18-30_AGO1_2

../src/flaimapper.py \
                     -v \
                     -o ../output/SRP006788/01.d_output_flaimapper.txt \
                     -f 1 \
                        ../share/small_RNA-seq_alignments/SRP006788/SRR207114_HeLa18-30_AGO1_2_RRP40

../src/flaimapper.py \
                     -v \
                     -o ../output/SRP006788/01.e_output_flaimapper.txt \
                     -f 1 \
                        ../share/small_RNA-seq_alignments/SRP006788/SRR207115_HeLa18-30_XRN1_2

../src/flaimapper.py \
                     -v \
                     -o ../output/SRP006788/01.f_output_flaimapper.txt \
                     -f 1 \
                        ../share/small_RNA-seq_alignments/SRP006788/SRR207116_HeLa18-30_N


### Step 02: aggregate data
grep -v Precursor ../output/SRP002175/01_output_flaimapper.txt > ../output/SRP002175/02_output_flaimapper.txt
cut -f 2 ../output/SRP002175/02_output_flaimapper.txt | sort | uniq > ../output/SRP002175/03_unique_ncRNAs_with_fragments.txt
grep -E "^>" ../share/annotations/ncRNA_annotation/ncRNdb09_with_tRNAs_and_Pseudogenes__21_oct_2011__hg19.fasta > ../output/SRP002175/04_unique_ncRNAs.txt
printf "[Fragments]\nMIRNA\t"$(grep =MIR ../output/SRP002175/02_output_flaimapper.txt | wc -l)"\nSNORD\t"$(grep SNORD ../output/SRP002175/02_output_flaimapper.txt | wc -l)"\nSNORA\t"$(grep SNORA ../output/SRP002175/02_output_flaimapper.txt | wc -l)"\nTRNA\t"$(grep NAME=TRNA ../output/SRP002175/02_output_flaimapper.txt | wc -l)"\nSCARNA\t"$(grep SCARNA ../output/SRP002175/02_output_flaimapper.txt | wc -l)"\nMISC\t"$(grep -v =MIR ../output/SRP002175/02_output_flaimapper.txt | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA  | wc -l)"\n\n[Affected ncRNAs/total reference ncRNAs]\nMIRNA\t"$(grep =MIR ../output/SRP002175/03_unique_ncRNAs_with_fragments.txt | wc -l)"/"$(grep =MIR      ../output/SRP002175/04_unique_ncRNAs.txt | wc -l)"\nSNORD\t"$(grep SNORD ../output/SRP002175/03_unique_ncRNAs_with_fragments.txt | wc -l)"/"$(grep SNORD     ../output/SRP002175/04_unique_ncRNAs.txt | wc -l)"\nSNORA\t"$(grep SNORA ../output/SRP002175/03_unique_ncRNAs_with_fragments.txt | wc -l)"/"$(grep SNORA     ../output/SRP002175/04_unique_ncRNAs.txt | wc -l)"\nTRNA\t"$(grep NAME=TRNA ../output/SRP002175/03_unique_ncRNAs_with_fragments.txt | wc -l)"/"$(grep NAME=TRNA ../output/SRP002175/04_unique_ncRNAs.txt | wc -l)"\nSCARNA\t"$(grep SCARNA ../output/SRP002175/03_unique_ncRNAs_with_fragments.txt | wc -l)"/"$(grep SCARNA    ../output/SRP002175/04_unique_ncRNAs.txt | wc -l)"\nMISC\t"$(grep -v =MIR ../output/SRP002175/03_unique_ncRNAs_with_fragments.txt | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA  | wc -l)"/"$(grep -v =MIR ../output/SRP002175/04_unique_ncRNAs.txt | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA  | wc -l)"\n\n" > ../output/SRP002175/05_summary.txt



SUFFIXES=( "a" "b" "c" "d" "e" "f" )
for SUFFIX in "${SUFFIXES[@]}"
do :
	grep -v Precursor "../output/SRP006788/01."$SUFFIX"_output_flaimapper.txt" > "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt"
	cut -f 2 "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt" | sort | uniq > "../output/SRP006788/03."$SUFFIX"_unique_ncRNAs_with_fragments.txt"
	grep -E "^>" ../share/annotations/ncRNA_annotation/ncRNdb09_with_tRNAs_and_Pseudogenes__21_oct_2011__hg19.fasta > ../output/SRP006788/04_unique_ncRNAs.txt
	printf "[Fragments]\nMIRNA\t"$(grep =MIR "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt" | wc -l)"\nSNORD\t"$(grep SNORD "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt" | wc -l)"\nSNORA\t"$(grep SNORA "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt" | wc -l)"\nTRNA\t"$(grep NAME=TRNA "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt" | wc -l)"\nSCARNA\t"$(grep SCARNA "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt" | wc -l)"\nMISC\t"$(grep -v =MISC "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt" | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA  | wc -l)"\n\n[Affected ncRNAs/total reference ncRNAs]\nMIRNA\t"$(grep =MIR "../output/SRP006788/03."$SUFFIX"_unique_ncRNAs_with_fragments.txt" | wc -l)"/"$(grep =MIR      ../output/SRP006788/04_unique_ncRNAs.txt | wc -l)"\nSNORD\t"$(grep SNORD "../output/SRP006788/03."$SUFFIX"_unique_ncRNAs_with_fragments.txt" | wc -l)"/"$(grep SNORD     ../output/SRP006788/04_unique_ncRNAs.txt | wc -l)"\nSNORA\t"$(grep SNORA "../output/SRP006788/03."$SUFFIX"_unique_ncRNAs_with_fragments.txt" | wc -l)"/"$(grep SNORA     ../output/SRP006788/04_unique_ncRNAs.txt | wc -l)"\nTRNA\t"$(grep NAME=TRNA "../output/SRP006788/03."$SUFFIX"_unique_ncRNAs_with_fragments.txt" | wc -l)"/"$(grep NAME=TRNA ../output/SRP006788/04_unique_ncRNAs.txt | wc -l)"\nSCARNA\t"$(grep SCARNA "../output/SRP006788/03."$SUFFIX"_unique_ncRNAs_with_fragments.txt" | wc -l)"/"$(grep SCARNA    ../output/SRP006788/04_unique_ncRNAs.txt | wc -l)"\nMISC\t"$(grep -v =MIR "../output/SRP006788/03."$SUFFIX"_unique_ncRNAs_with_fragments.txt" | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA  | wc -l)"/"$(grep -v =MIR "../output/SRP002175/04_unique_ncRNAs.txt" | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA  | wc -l)"\n" > "../output/SRP006788/05."$SUFFIX"_summary.txt"
done



# Fig 1:
cd figures/fig_1
python alignment_plot.py > /dev/null
R --vanilla --slave -f "alignment_plot.R" > /dev/null
cd ../../

# Fig 4:
cd figures/fig_4
R --vanilla --slave -f "filter_shelf.R"
cd ../../

# Fig 5:
cd figures/fig_5

cd 0_top_validation_FlaiMapper_miRBase
python 0_top_validation_FlaiMapper_miRBase_SRP002175.py
R --vanilla --slave -f "0_top_validation_FlaiMapper_miRBase_SRP002175.R"

python validation_FlaiMapper_miRBase_SRP006788.py
R --vanilla --slave -f "validation_FlaiMapper_miRBase_SRP006788.R"
cd ..

cd 1_bottom_validation_BlockBuster_miRBase

echo "Running BlockBuster (will take some time)"
#./validation_BlockBuster_miRBase_SRP002175.sh

echo "Assessing BlockBuster (will also take some time)"
##python validation_BlockBuster_miRBase_SRP002175__optimal_settings.py
R --vanilla --slave -f "validation_BlockBuster_miRBase_SRP002175__optimal_settings.R"

python 1_bottom_validation_BlockBuster_miRBase_SRP002175.py
R --vanilla --slave -f "1_bottom_validation_BlockBuster_miRBase_SRP002175.R"
cd ..

cd ../../

## Fig 6:
cd figures/fig_6

python validation_miRBase_SRP002175__sequencing_depth_vs_offset.py
R --vanilla --slave -f "validation_miRBase_SRP002175__sequencing_depth_vs_offset.R"

cd ../../

# Fig 7:
cd figures/fig_7
R --vanilla --slave -f "PCA.R"
cd ../../

# Fig 8:
cd figures/fig_8
python extract_motifs.py
cd ../../