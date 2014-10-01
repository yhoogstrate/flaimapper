#!/bin/bash

echo "01. Unpacking alignments"

cd "../share/small_RNA-seq_alignments/"


# ---------------------------SRP002175----------------------------------

cd "SRP002175"

echo "  - Dataset: SRP002175"

SAMPLES=( "SRR038852" "SRR038853" "SRR038854" "SRR038855" "SRR038856" "SRR038857" "SRR038858" "SRR038859" "SRR038860" "SRR038861" "SRR038862" "SRR038863" )

for SAMPLE in ${SAMPLES[@]}
do :
	if [ -e $SAMPLE".tar.gz" ]
	then
		echo "    * Extracting: $SAMPLE"
		tar xzvf $SAMPLE".tar.gz" > /dev/null
	else
		echo "    ! Could not find archive: "$SAMPLE".tar.gz"
	fi
done

unset SAMPLE SAMPLES

cd ../

# ----------------------------------------------------------------------

# ---------------------------SRP006788----------------------------------

cd "SRP006788"

echo "  - Dataset: SRP006788"

SAMPLES=( "SRR207111_HeLa18-30" "SRR207112_HeLa18-30_RRP40" "SRR207113_HeLa18-30_AGO1_2" "SRR207114_HeLa18-30_AGO1_2_RRP40" "SRR207115_HeLa18-30_XRN1_2" "SRR207116_HeLa18-30_N" )

for SAMPLE in ${SAMPLES[@]}
do :
	if [ -e $SAMPLE".tar.gz" ]
	then
		echo "    * Extracting: $SAMPLE"
		tar xzvf $SAMPLE".tar.gz" > /dev/null
	else
		echo "    ! Could not find archive: "$SAMPLE".tar.gz"
	fi
done

unset SAMPLE SAMPLES

cd ../

# ----------------------------------------------------------------------

# ---------------------------SRP028959----------------------------------

cd "SRP028959"

echo "  - Dataset: SRP028959"

SAMPLES=( "SRR954957" "SRR954958" "SRR954959" )

for SAMPLE in ${SAMPLES[@]}
do :
	if [ -e $SAMPLE".tar.gz" ]
	then
		echo "    * Extracting: $SAMPLE"
		tar xzvf $SAMPLE".tar.gz" > /dev/null
	else
		echo "    ! Could not find archive: "$SAMPLE".tar.gz"
	fi
done

unset SAMPLE SAMPLES

cd ../

# ----------------------------------------------------------------------

# ---------------------------SRP034013----------------------------------

cd "SRP034013"

echo "  - Dataset: SRP034013"

SAMPLES=( "SRR1049397" "SRR1049398" "SRR1049399" )

for SAMPLE in ${SAMPLES[@]}
do :
	if [ -e $SAMPLE".tar.gz" ]
	then
		echo "    * Extracting: $SAMPLE"
		tar xzvf $SAMPLE".tar.gz" > /dev/null
	else
		echo "    ! Could not find archive: "$SAMPLE".tar.gz"
	fi
done

unset SAMPLE SAMPLES

cd ../

# ----------------------------------------------------------------------

# ---------------------------SRP041082----------------------------------

cd "SRP041082"

echo "  - Dataset: SRP041082"

SAMPLES=( "SRR1232072" "SRR1232073" )

for SAMPLE in ${SAMPLES[@]}
do :
	if [ -e $SAMPLE".tar.gz" ]
	then
		echo "    * Extracting: $SAMPLE"
		tar xzvf $SAMPLE".tar.gz" > /dev/null
	else
		echo "    ! Could not find archive: "$SAMPLE".tar.gz"
	fi
done

unset SAMPLE SAMPLES

cd ../

# ----------------------------------------------------------------------


cd ../../scripts
