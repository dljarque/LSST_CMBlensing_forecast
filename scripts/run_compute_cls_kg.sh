#!/bin/bash

PYTHON=/home/dljarque/miniforge3/envs/IFAE/bin/python
SCRIPT=./scripts/compute_colore_cls.py

$PYTHON $SCRIPT 'cmblensing_x_sources' \
\
--kappa_map ./data/test_box_4/results/out_kappa_z000.fits \
--sources_prefix ./data/test_box_4/results/out_srcs_s5 \
\
--lmin 2 \
--nbins 30 \
--out_dir ./data/test_box_4/cls \
--filename 'desi_qso' \
#--zmean 2.0 \
#--dz 0.2 \
#--sigma_z 0.05