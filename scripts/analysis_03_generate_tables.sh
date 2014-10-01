#!/bin/bash

echo "03. Generate tables"

# ---------------------------SRP002175----------------------------------

echo "  - Dataset: SRP002175 (experiments merged)"

FILE1="../output/flaimapper/SRP002175_merged/01_output_flaimapper.txt"
FILE2="../output/flaimapper/SRP002175_merged/02_output_flaimapper.txt"
FILE3="../output/flaimapper/SRP002175_merged/03_unique_ncRNAs_with_fragments.txt"
FILE4="../output/flaimapper/SRP002175_merged/04_unique_ncRNAs.txt"
FILE5="../output/flaimapper/SRP002175_merged/05_summary.txt"

grep -v Precursor $FILE1 > $FILE2
cut -f 2 $FILE2 | sort | uniq >  $FILE3
grep -E "^>" ../share/annotations/ncRNA_annotation/ncrnadb09.fa > $FILE4

printf "[Fragments]\n" > $FILE5
printf  "MIRNA\t"$(grep =MIR $FILE2      | wc -l)"\n" >> $FILE5
printf  "SNORD\t"$(grep SNORD $FILE2     | wc -l)"\n" >> $FILE5
printf  "SNORA\t"$(grep SNORA $FILE2     | wc -l)"\n" >> $FILE5
printf   "TRNA\t"$(grep NAME=TRNA $FILE2 | wc -l)"\n" >> $FILE5
printf "SCARNA\t"$(grep SCARNA $FILE2    | wc -l)"\n" >> $FILE5
printf   "MISC\t"$(grep -v =MIR $FILE2   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA  | wc -l)"\n" >> $FILE5
printf "\n" >> $FILE5
printf "[Affected ncRNAs/total reference ncRNAs]\n" >> $FILE5
printf  "MIRNA\t"$(grep =MIR $FILE3      | wc -l)"/"$(grep =MIR $FILE4 | wc -l)"\n" >> $FILE5
printf  "SNORD\t"$(grep SNORD $FILE3     | wc -l)"/"$(grep SNORD $FILE4 | wc -l)"\n" >> $FILE5
printf  "SNORA\t"$(grep SNORA $FILE3     | wc -l)"/"$(grep SNORA  $FILE4 | wc -l)"\n" >> $FILE5
printf   "TRNA\t"$(grep NAME=TRNA $FILE3 | wc -l)"/"$(grep NAME=TRNA $FILE4 | wc -l)"\n" >> $FILE5
printf "SCARNA\t"$(grep SCARNA $FILE3    | wc -l)"/"$(grep SCARNA $FILE4 | wc -l)"\n" >> $FILE5
printf   "MISC\t"$(grep -v =MIR $FILE3   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA  | wc -l)"/"$(grep -v =MIR $FILE4 | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA  | wc -l)"\n" >> $FILE5
printf "\n" >> $FILE5

# ----------------------------------------------------------------------

echo "  - Dataset: SRP002175"

SAMPLES=( "SRR038852" "SRR038853" "SRR038854" "SRR038855" "SRR038856" "SRR038857" "SRR038858" "SRR038859" "SRR038860" "SRR038861" "SRR038862" "SRR038863" )

for SAMPLE in ${SAMPLES[@]}
do :
	echo "    * Experiment "$SAMPLE
	
	FILE1="../output/flaimapper/SRP002175/01_output_flaimapper_"$SAMPLE".txt"
	FILE2="../output/flaimapper/SRP002175/02_output_flaimapper_"$SAMPLE".txt"
	FILE3="../output/flaimapper/SRP002175/03_unique_ncRNAs_with_fragments_"$SAMPLE".txt"
	FILE4="../output/flaimapper/SRP002175/04_unique_ncRNAs_"$SAMPLE".txt"
	FILE5="../output/flaimapper/SRP002175/05_summary_"$SAMPLE".txt"
	
	grep -v Precursor $FILE1 > $FILE2
	cut -f 2 $FILE2 | sort | uniq >  $FILE3
	grep -E "^>" ../share/annotations/ncRNA_annotation/ncrnadb09.fa > $FILE4
	
	printf "[Fragments]\n" > $FILE5
	printf  "MIRNA\t"$(grep =MIR $FILE2      | wc -l)"\n" >> $FILE5
	printf  "SNORD\t"$(grep SNORD $FILE2     | wc -l)"\n" >> $FILE5
	printf  "SNORA\t"$(grep SNORA $FILE2     | wc -l)"\n" >> $FILE5
	printf   "TRNA\t"$(grep NAME=TRNA $FILE2 | wc -l)"\n" >> $FILE5
	printf "SCARNA\t"$(grep SCARNA $FILE2    | wc -l)"\n" >> $FILE5
	printf   "MISC\t"$(grep -v =MIR $FILE2   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA  | wc -l)"\n" >> $FILE5
	printf "\n" >> $FILE5
	printf "[Affected ncRNAs/total reference ncRNAs]\n" >> $FILE5
	printf  "MIRNA\t"$(grep =MIR $FILE3      | wc -l)"/"$(grep =MIR $FILE4 | wc -l)"\n" >> $FILE5
	printf  "SNORD\t"$(grep SNORD $FILE3     | wc -l)"/"$(grep SNORD $FILE4 | wc -l)"\n" >> $FILE5
	printf  "SNORA\t"$(grep SNORA $FILE3     | wc -l)"/"$(grep SNORA  $FILE4 | wc -l)"\n" >> $FILE5
	printf   "TRNA\t"$(grep NAME=TRNA $FILE3 | wc -l)"/"$(grep NAME=TRNA $FILE4 | wc -l)"\n" >> $FILE5
	printf "SCARNA\t"$(grep SCARNA $FILE3    | wc -l)"/"$(grep SCARNA $FILE4 | wc -l)"\n" >> $FILE5
	printf   "MISC\t"$(grep -v =MIR $FILE3   | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA  | wc -l)"/"$(grep -v =MIR $FILE4 | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA  | wc -l)"\n" >> $FILE5
	printf "\n" >> $FILE5
done

unset SAMPLE SAMPLES

# ----------------------------------------------------------------------

# ---------------------------SRP006788----------------------------------

echo "  - Dataset: SRP006788"

#SAMPLES=( "a" "b" "c" "d" "e" "f" )
#for SAMPLE in ${SAMPLES[@]}
#do :
	#grep -v Precursor "../output/SRP006788/01."$SUFFIX"_output_flaimapper.txt" > "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt"
	#cut -f 2 "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt" | sort | uniq > "../output/SRP006788/03."$SUFFIX"_unique_ncRNAs_with_fragments.txt"
	#grep -E "^>" ../share/annotations/ncRNA_annotation/ncrnadb09.fa > ../output/SRP006788/04_unique_ncRNAs.txt
	#printf "[Fragments]\nMIRNA\t"$(grep =MIR "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt" | wc -l)"\nSNORD\t"$(grep SNORD "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt" | wc -l)"\nSNORA\t"$(grep SNORA "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt" | wc -l)"\nTRNA\t"$(grep NAME=TRNA "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt" | wc -l)"\nSCARNA\t"$(grep SCARNA "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt" | wc -l)"\nMISC\t"$(grep -v =MISC "../output/SRP006788/02."$SUFFIX"_output_flaimapper.txt" | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA  | wc -l)"\n\n[Affected ncRNAs/total reference ncRNAs]\nMIRNA\t"$(grep =MIR "../output/SRP006788/03."$SUFFIX"_unique_ncRNAs_with_fragments.txt" | wc -l)"/"$(grep =MIR      ../output/SRP006788/04_unique_ncRNAs.txt | wc -l)"\nSNORD\t"$(grep SNORD "../output/SRP006788/03."$SUFFIX"_unique_ncRNAs_with_fragments.txt" | wc -l)"/"$(grep SNORD     ../output/SRP006788/04_unique_ncRNAs.txt | wc -l)"\nSNORA\t"$(grep SNORA "../output/SRP006788/03."$SUFFIX"_unique_ncRNAs_with_fragments.txt" | wc -l)"/"$(grep SNORA     ../output/SRP006788/04_unique_ncRNAs.txt | wc -l)"\nTRNA\t"$(grep NAME=TRNA "../output/SRP006788/03."$SUFFIX"_unique_ncRNAs_with_fragments.txt" | wc -l)"/"$(grep NAME=TRNA ../output/SRP006788/04_unique_ncRNAs.txt | wc -l)"\nSCARNA\t"$(grep SCARNA "../output/SRP006788/03."$SUFFIX"_unique_ncRNAs_with_fragments.txt" | wc -l)"/"$(grep SCARNA    ../output/SRP006788/04_unique_ncRNAs.txt | wc -l)"\nMISC\t"$(grep -v =MIR "../output/SRP006788/03."$SUFFIX"_unique_ncRNAs_with_fragments.txt" | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA  | wc -l)"/"$(grep -v =MIR "../output/SRP002175/04_unique_ncRNAs.txt" | grep -v SNORD | grep -v SNORA | grep -v NAME=TRNA | grep -v SCARNA  | wc -l)"\n" > "../output/SRP006788/05."$SUFFIX"_summary.txt"
#done

# ----------------------------------------------------------------------

# ---------------------------SRP028959----------------------------------

echo "  - Dataset: SRP028959"


# ----------------------------------------------------------------------

# ---------------------------SRP034013----------------------------------

echo "  - Dataset: SRP034013 (merged experiments)"



# ----------------------------------------------------------------------

# ---------------------------SRP041082----------------------------------

echo "  - Dataset: SRP041082 (merged experiments)"


# ----------------------------------------------------------------------
