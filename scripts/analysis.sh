#!/bin/bash

source analysis_01_extract_alignments.sh
source analysis_02_run_flaimapper.sh
source analysis_03_generate_tables.sh


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
