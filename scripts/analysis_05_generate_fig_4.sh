#!/bin/bash

# ----------------------------------------------------------------------
echo "05. Generate Fig. 4"
cd figures/fig_4
R --vanilla --slave -f "filter_shelf.R" '--args /tmp ../../../output/figures/fig_4/' > /dev/null
cd ../../
# ----------------------------------------------------------------------
