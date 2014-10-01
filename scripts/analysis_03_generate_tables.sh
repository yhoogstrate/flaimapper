#!/bin/bash

echo "03. Generate tables"

# ---------------------------SRP002175----------------------------------

echo "  - Dataset: SRP002175 (experiments merged)"

ROOT_DIR="../output/flaimapper/SRP002175_merged"

FILE1=$ROOT_DIR"/01_output_flaimapper.txt"
FILE2=$ROOT_DIR"/02_output_flaimapper.txt"
FILE3=$ROOT_DIR"/03_unique_ncRNAs_with_fragments.txt"
FILE4=$ROOT_DIR"/04_unique_ncRNAs.txt"
FILE5=$ROOT_DIR"/05_summary.txt"

grep -v Precursor $FILE1 > $FILE2
cut -f 2 $FILE2 | sort | uniq >  $FILE3
grep -E "^>" ../share/annotations/ncRNA_annotation/ncrnadb09.fa > $FILE4

printf "[Fragments]\n" > $FILE5
printf  "MIRNA\t"$(grep =MIR $FILE2      | wc -l)"\n" >> $FILE5
printf  "SNORD\t"$(grep SNORD $FILE2     | wc -l)"\n" >> $FILE5
printf  "SNORA\t"$(grep SNORA $FILE2     | wc -l)"\n" >> $FILE5
printf   "TRNA\t"$(grep NAME=TRNA $FILE2 | wc -l)"\n" >> $FILE5
printf "SCARNA\t"$(grep SCARNA $FILE2    | wc -l)"\n" >> $FILE5
printf   "MISC\t"$(grep -v =MIR $FILE2   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"\n" >> $FILE5
printf "\n" >> $FILE5
printf "[Affected ncRNAs/total reference ncRNAs]\n" >> $FILE5
printf  "MIRNA\t"$(grep =MIR $FILE3      | wc -l)"/"$(grep =MIR $FILE4 | wc -l)"\n" >> $FILE5
printf  "SNORD\t"$(grep SNORD $FILE3     | wc -l)"/"$(grep SNORD $FILE4 | wc -l)"\n" >> $FILE5
printf  "SNORA\t"$(grep SNORA $FILE3     | wc -l)"/"$(grep SNORA  $FILE4 | wc -l)"\n" >> $FILE5
printf   "TRNA\t"$(grep NAME=TRNA $FILE3 | wc -l)"/"$(grep NAME=TRNA $FILE4 | wc -l)"\n" >> $FILE5
printf "SCARNA\t"$(grep SCARNA $FILE3    | wc -l)"/"$(grep SCARNA $FILE4 | wc -l)"\n" >> $FILE5
printf   "MISC\t"$(grep -v =MIR $FILE3   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"/"$(grep -v =MIR $FILE4 | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"\n" >> $FILE5
printf "\n" >> $FILE5

unset ROOT_DIR FILE1 FILE2 FILE3 FILE4 FILE5

# ----------------------------------------------------------------------

echo "  - Dataset: SRP002175"

ROOT_DIR="../output/flaimapper/SRP002175"
SAMPLES=( "SRR038852" "SRR038853" "SRR038854" "SRR038855" "SRR038856" "SRR038857" "SRR038858" "SRR038859" "SRR038860" "SRR038861" "SRR038862" "SRR038863" )

for SAMPLE in ${SAMPLES[@]}
do :
	echo "    * Experiment "$SAMPLE
	
	FILE1=$ROOT_DIR"/01_output_flaimapper_"$SAMPLE".txt"
	FILE2=$ROOT_DIR"/02_output_flaimapper_"$SAMPLE".txt"
	FILE3=$ROOT_DIR"/03_unique_ncRNAs_with_fragments_"$SAMPLE".txt"
	FILE4=$ROOT_DIR"/04_unique_ncRNAs_"$SAMPLE".txt"
	FILE5=$ROOT_DIR"/05_summary_"$SAMPLE".txt"
	
	grep -v Precursor $FILE1 > $FILE2
	cut -f 2 $FILE2 | sort | uniq >  $FILE3
	grep -E "^>" ../share/annotations/ncRNA_annotation/ncrnadb09.fa > $FILE4
	
	printf "[Fragments]\n" > $FILE5
	printf  "MIRNA\t"$(grep =MIR $FILE2      | wc -l)"\n" >> $FILE5
	printf  "SNORD\t"$(grep SNORD $FILE2     | wc -l)"\n" >> $FILE5
	printf  "SNORA\t"$(grep SNORA $FILE2     | wc -l)"\n" >> $FILE5
	printf   "TRNA\t"$(grep NAME=TRNA $FILE2 | wc -l)"\n" >> $FILE5
	printf "SCARNA\t"$(grep SCARNA $FILE2    | wc -l)"\n" >> $FILE5
	printf   "MISC\t"$(grep -v =MIR $FILE2   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"\n" >> $FILE5
	printf "\n" >> $FILE5
	printf "[Affected ncRNAs/total reference ncRNAs]\n" >> $FILE5
	printf  "MIRNA\t"$(grep =MIR $FILE3      | wc -l)"/"$(grep =MIR $FILE4 | wc -l)"\n" >> $FILE5
	printf  "SNORD\t"$(grep SNORD $FILE3     | wc -l)"/"$(grep SNORD $FILE4 | wc -l)"\n" >> $FILE5
	printf  "SNORA\t"$(grep SNORA $FILE3     | wc -l)"/"$(grep SNORA  $FILE4 | wc -l)"\n" >> $FILE5
	printf   "TRNA\t"$(grep NAME=TRNA $FILE3 | wc -l)"/"$(grep NAME=TRNA $FILE4 | wc -l)"\n" >> $FILE5
	printf "SCARNA\t"$(grep SCARNA $FILE3    | wc -l)"/"$(grep SCARNA $FILE4 | wc -l)"\n" >> $FILE5
	printf   "MISC\t"$(grep -v =MIR $FILE3   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"/"$(grep -v =MIR $FILE4 | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"\n" >> $FILE5
	printf "\n" >> $FILE5
done

unset SAMPLE SAMPLES ROOT_DIR FILE1 FILE2 FILE3 FILE4 FILE5

# ----------------------------------------------------------------------

# ---------------------------SRP006788----------------------------------

echo "  - Dataset: SRP006788"

ROOT_DIR="../output/flaimapper/SRP006788"
SAMPLES=( "SRR207111_HeLa18-30" "SRR207112_HeLa18-30_RRP40" "SRR207113_HeLa18-30_AGO1_2" "SRR207114_HeLa18-30_AGO1_2_RRP40" "SRR207115_HeLa18-30_XRN1_2" "SRR207116_HeLa18-30_N" )

for SAMPLE in ${SAMPLES[@]}
do :
	echo "    * Experiment "$SAMPLE
	
	FILE1=$ROOT_DIR"/01_output_flaimapper_"$SAMPLE".txt"
	FILE2=$ROOT_DIR"/02_output_flaimapper_"$SAMPLE".txt"
	FILE3=$ROOT_DIR"/03_unique_ncRNAs_with_fragments_"$SAMPLE".txt"
	FILE4=$ROOT_DIR"/04_unique_ncRNAs_"$SAMPLE".txt"
	FILE5=$ROOT_DIR"/05_summary_"$SAMPLE".txt"
	
	grep -v Precursor $FILE1 > $FILE2
	cut -f 2 $FILE2 | sort | uniq >  $FILE3
	grep -E "^>" ../share/annotations/ncRNA_annotation/ncrnadb09.fa > $FILE4
	
	printf "[Fragments]\n" > $FILE5
	printf  "MIRNA\t"$(grep =MIR $FILE2      | wc -l)"\n" >> $FILE5
	printf  "SNORD\t"$(grep SNORD $FILE2     | wc -l)"\n" >> $FILE5
	printf  "SNORA\t"$(grep SNORA $FILE2     | wc -l)"\n" >> $FILE5
	printf   "TRNA\t"$(grep NAME=TRNA $FILE2 | wc -l)"\n" >> $FILE5
	printf "SCARNA\t"$(grep SCARNA $FILE2    | wc -l)"\n" >> $FILE5
	printf   "MISC\t"$(grep -v =MIR $FILE2   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"\n" >> $FILE5
	printf "\n" >> $FILE5
	printf "[Affected ncRNAs/total reference ncRNAs]\n" >> $FILE5
	printf  "MIRNA\t"$(grep =MIR $FILE3      | wc -l)"/"$(grep =MIR $FILE4 | wc -l)"\n" >> $FILE5
	printf  "SNORD\t"$(grep SNORD $FILE3     | wc -l)"/"$(grep SNORD $FILE4 | wc -l)"\n" >> $FILE5
	printf  "SNORA\t"$(grep SNORA $FILE3     | wc -l)"/"$(grep SNORA  $FILE4 | wc -l)"\n" >> $FILE5
	printf   "TRNA\t"$(grep NAME=TRNA $FILE3 | wc -l)"/"$(grep NAME=TRNA $FILE4 | wc -l)"\n" >> $FILE5
	printf "SCARNA\t"$(grep SCARNA $FILE3    | wc -l)"/"$(grep SCARNA $FILE4 | wc -l)"\n" >> $FILE5
	printf   "MISC\t"$(grep -v =MIR $FILE3   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"/"$(grep -v =MIR $FILE4 | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"\n" >> $FILE5
	printf "\n" >> $FILE5
done

unset SAMPLE SAMPLES FILE1 FILE2 FILE3 FILE4 FILE5

# ----------------------------------------------------------------------

# ---------------------------SRP028959----------------------------------

echo "  - Dataset: SRP028959"

ROOT_DIR="../output/flaimapper/SRP028959"
SAMPLES=( "SRR954957" "SRR954958" "SRR954959" )

for SAMPLE in ${SAMPLES[@]}
do :
	echo "    * Experiment "$SAMPLE
	
	FILE1=$ROOT_DIR"/01_output_flaimapper_"$SAMPLE".txt"
	FILE2=$ROOT_DIR"/02_output_flaimapper_"$SAMPLE".txt"
	FILE3=$ROOT_DIR"/03_unique_ncRNAs_with_fragments_"$SAMPLE".txt"
	FILE4=$ROOT_DIR"/04_unique_ncRNAs_"$SAMPLE".txt"
	FILE5=$ROOT_DIR"/05_summary_"$SAMPLE".txt"
	
	grep -v Precursor $FILE1 > $FILE2
	cut -f 2 $FILE2 | sort | uniq >  $FILE3
	grep -E "^>" ../share/annotations/ncRNA_annotation/ncrnadb09.fa > $FILE4
	
	printf "[Fragments]\n" > $FILE5
	printf  "MIRNA\t"$(grep =MIR $FILE2      | wc -l)"\n" >> $FILE5
	printf  "SNORD\t"$(grep SNORD $FILE2     | wc -l)"\n" >> $FILE5
	printf  "SNORA\t"$(grep SNORA $FILE2     | wc -l)"\n" >> $FILE5
	printf   "TRNA\t"$(grep NAME=TRNA $FILE2 | wc -l)"\n" >> $FILE5
	printf "SCARNA\t"$(grep SCARNA $FILE2    | wc -l)"\n" >> $FILE5
	printf   "MISC\t"$(grep -v =MIR $FILE2   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"\n" >> $FILE5
	printf "\n" >> $FILE5
	printf "[Affected ncRNAs/total reference ncRNAs]\n" >> $FILE5
	printf  "MIRNA\t"$(grep =MIR $FILE3      | wc -l)"/"$(grep =MIR $FILE4 | wc -l)"\n" >> $FILE5
	printf  "SNORD\t"$(grep SNORD $FILE3     | wc -l)"/"$(grep SNORD $FILE4 | wc -l)"\n" >> $FILE5
	printf  "SNORA\t"$(grep SNORA $FILE3     | wc -l)"/"$(grep SNORA  $FILE4 | wc -l)"\n" >> $FILE5
	printf   "TRNA\t"$(grep NAME=TRNA $FILE3 | wc -l)"/"$(grep NAME=TRNA $FILE4 | wc -l)"\n" >> $FILE5
	printf "SCARNA\t"$(grep SCARNA $FILE3    | wc -l)"/"$(grep SCARNA $FILE4 | wc -l)"\n" >> $FILE5
	printf   "MISC\t"$(grep -v =MIR $FILE3   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"/"$(grep -v =MIR $FILE4 | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"\n" >> $FILE5
	printf "\n" >> $FILE5
done

unset SAMPLE SAMPLES FILE1 FILE2 FILE3 FILE4 FILE5

# ----------------------------------------------------------------------

# ---------------------------SRP034013----------------------------------

echo "  - Dataset: SRP034013 (merged experiments)"

ROOT_DIR="../output/flaimapper/SRP034013_merged"

FILE1=$ROOT_DIR"/01_output_flaimapper.txt"
FILE2=$ROOT_DIR"/02_output_flaimapper.txt"
FILE3=$ROOT_DIR"/03_unique_ncRNAs_with_fragments.txt"
FILE4=$ROOT_DIR"/04_unique_ncRNAs.txt"
FILE5=$ROOT_DIR"/05_summary.txt"

grep -v Precursor $FILE1 > $FILE2
cut -f 2 $FILE2 | sort | uniq >  $FILE3
grep -E "^>" ../share/annotations/ncRNA_annotation/ncrnadb09.fa > $FILE4

printf "[Fragments]\n" > $FILE5
printf  "MIRNA\t"$(grep =MIR $FILE2      | wc -l)"\n" >> $FILE5
printf  "SNORD\t"$(grep SNORD $FILE2     | wc -l)"\n" >> $FILE5
printf  "SNORA\t"$(grep SNORA $FILE2     | wc -l)"\n" >> $FILE5
printf   "TRNA\t"$(grep NAME=TRNA $FILE2 | wc -l)"\n" >> $FILE5
printf "SCARNA\t"$(grep SCARNA $FILE2    | wc -l)"\n" >> $FILE5
printf   "MISC\t"$(grep -v =MIR $FILE2   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"\n" >> $FILE5
printf "\n" >> $FILE5
printf "[Affected ncRNAs/total reference ncRNAs]\n" >> $FILE5
printf  "MIRNA\t"$(grep =MIR $FILE3      | wc -l)"/"$(grep =MIR $FILE4 | wc -l)"\n" >> $FILE5
printf  "SNORD\t"$(grep SNORD $FILE3     | wc -l)"/"$(grep SNORD $FILE4 | wc -l)"\n" >> $FILE5
printf  "SNORA\t"$(grep SNORA $FILE3     | wc -l)"/"$(grep SNORA  $FILE4 | wc -l)"\n" >> $FILE5
printf   "TRNA\t"$(grep NAME=TRNA $FILE3 | wc -l)"/"$(grep NAME=TRNA $FILE4 | wc -l)"\n" >> $FILE5
printf "SCARNA\t"$(grep SCARNA $FILE3    | wc -l)"/"$(grep SCARNA $FILE4 | wc -l)"\n" >> $FILE5
printf   "MISC\t"$(grep -v =MIR $FILE3   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"/"$(grep -v =MIR $FILE4 | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"\n" >> $FILE5
printf "\n" >> $FILE5

unset ROOT_DIR FILE1 FILE2 FILE3 FILE4 FILE5

# ----------------------------------------------------------------------

# ---------------------------SRP041082----------------------------------

echo "  - Dataset: SRP041082 (merged experiments)"

ROOT_DIR="../output/flaimapper/SRP041082_merged"

FILE1=$ROOT_DIR"/01_output_flaimapper.txt"
FILE2=$ROOT_DIR"/02_output_flaimapper.txt"
FILE3=$ROOT_DIR"/03_unique_ncRNAs_with_fragments.txt"
FILE4=$ROOT_DIR"/04_unique_ncRNAs.txt"
FILE5=$ROOT_DIR"/05_summary.txt"

grep -v Precursor $FILE1 > $FILE2
cut -f 2 $FILE2 | sort | uniq >  $FILE3
grep -E "^>" ../share/annotations/ncRNA_annotation/ncrnadb09.fa > $FILE4

printf "[Fragments]\n" > $FILE5
printf  "MIRNA\t"$(grep =MIR $FILE2      | wc -l)"\n" >> $FILE5
printf  "SNORD\t"$(grep SNORD $FILE2     | wc -l)"\n" >> $FILE5
printf  "SNORA\t"$(grep SNORA $FILE2     | wc -l)"\n" >> $FILE5
printf   "TRNA\t"$(grep NAME=TRNA $FILE2 | wc -l)"\n" >> $FILE5
printf "SCARNA\t"$(grep SCARNA $FILE2    | wc -l)"\n" >> $FILE5
printf   "MISC\t"$(grep -v =MIR $FILE2   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"\n" >> $FILE5
printf "\n" >> $FILE5
printf "[Affected ncRNAs/total reference ncRNAs]\n" >> $FILE5
printf  "MIRNA\t"$(grep =MIR $FILE3      | wc -l)"/"$(grep =MIR $FILE4 | wc -l)"\n" >> $FILE5
printf  "SNORD\t"$(grep SNORD $FILE3     | wc -l)"/"$(grep SNORD $FILE4 | wc -l)"\n" >> $FILE5
printf  "SNORA\t"$(grep SNORA $FILE3     | wc -l)"/"$(grep SNORA  $FILE4 | wc -l)"\n" >> $FILE5
printf   "TRNA\t"$(grep NAME=TRNA $FILE3 | wc -l)"/"$(grep NAME=TRNA $FILE4 | wc -l)"\n" >> $FILE5
printf "SCARNA\t"$(grep SCARNA $FILE3    | wc -l)"/"$(grep SCARNA $FILE4 | wc -l)"\n" >> $FILE5
printf   "MISC\t"$(grep -v =MIR $FILE3   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"/"$(grep -v =MIR $FILE4 | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA | wc -l)"\n" >> $FILE5
printf "\n" >> $FILE5

unset ROOT_DIR FILE1 FILE2 FILE3 FILE4 FILE5

# ----------------------------------------------------------------------
