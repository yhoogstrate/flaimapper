#!/bin/bash

# ----------------------------------------------------------------------
echo "04. Generate Fig. 1"
cd figures/fig_1
python alignment_plot.py > /dev/null
R --vanilla --slave -f "alignment_plot.R" '--args /tmp ../../../output-sslm/figures/fig_1/' > /dev/null
cd ../../
# ----------------------------------------------------------------------
