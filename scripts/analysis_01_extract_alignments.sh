#!/bin/bash

echo "01. Unpacking alignments"

cd "../share/small_RNA-seq_alignments/"


# ---------------------------SRP002175----------------------------------

cd "SRP002175"

echo "  - Dataset: SRP002175"

FILES=( "SRR038852.tar.gz" "SRR038853.tar.gz" "SRR038854.tar.gz" "SRR038855.tar.gz" "SRR038856.tar.gz" "SRR038857.tar.gz" "SRR038858.tar.gz" "SRR038859.tar.gz" "SRR038860.tar.gz" "SRR038861.tar.gz" "SRR038862.tar.gz" "SRR038863.tar.gz" )

for FILE in "${FILES[@]}"
do :

  if [ -e $FILE ]
  then
    echo "    * Extracting: $FILE"
    tar xzvf $FILE > /dev/null
  else
    echo "    ! Could not find archive: $FILE"
  fi

done

unset FILE FILES

cd ../

# ----------------------------------------------------------------------

# ---------------------------SRP006788----------------------------------

cd "SRP006788"

echo "  - Dataset: SRP006788"

FILES=( "SRR207111_HeLa18-30.tar.gz" "SRR207112_HeLa18-30_RRP40.tar.gz" "SRR207113_HeLa18-30_AGO1_2.tar.gz" "SRR207114_HeLa18-30_AGO1_2_RRP40.tar.gz" "SRR207115_HeLa18-30_XRN1_2.tar.gz" "SRR207116_HeLa18-30_N.tar.gz" )

for FILE in "${FILES[@]}"
do :

  if [ -e $FILE ]
  then
    echo "    * Extracting: $FILE"
    tar xzvf $FILE > /dev/null
  else
    echo "    ! Could not find archive: $FILE"
  fi

done

unset FILE FILES

cd ../

# ----------------------------------------------------------------------

# ---------------------------SRP028959----------------------------------

cd "SRP028959"

echo "  - Dataset: SRP028959"

FILES=( "SRR954957.tar.gz" "SRR954958.tar.gz" "SRR954959.tar.gz" )

for FILE in "${FILES[@]}"
do :

  if [ -e $FILE ]
  then
    echo "    * Extracting: $FILE"
    tar xzvf $FILE > /dev/null
  else
    echo "    ! Could not find archive: $FILE"
  fi

done

unset FILE FILES

cd ../

# ----------------------------------------------------------------------

# ---------------------------SRP034013----------------------------------

cd "SRP034013"

echo "  - Dataset: SRP034013"

FILES=( "SRR1049397.tar.gz" "SRR1049398.tar.gz" "SRR1049399.tar.gz" )

for FILE in "${FILES[@]}"
do :

  if [ -e $FILE ]
  then
    echo "    * Extracting: $FILE"
    tar xzvf $FILE > /dev/null
  else
    echo "    ! Could not find archive: $FILE"
  fi

done

unset FILE FILES

cd ../

# ----------------------------------------------------------------------

# ---------------------------SRP041082----------------------------------

cd "SRP041082"

echo "  - Dataset: SRP041082"

FILES=( "SRR1232072.tar.gz" "SRR1232073.tar.gz" )

for FILE in "${FILES[@]}"
do :

  if [ -e $FILE ]
  then
    echo "    * Extracting: $FILE"
    tar xzvf $FILE > /dev/null
  else
    echo "    ! Could not find archive: $FILE"
  fi

done

unset FILE FILES

cd ../

# ----------------------------------------------------------------------


cd ../../scripts
