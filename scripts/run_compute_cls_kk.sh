#!/bin/bash

PYTHON=/home/dljarque/miniforge3/envs/IFAE/bin/python
SCRIPT=./scripts/compute_colore_cls.py

$PYTHON $SCRIPT 'cmb_lensing_auto' \
\
--kappa_map ./data/test_box_4/results/out_kappa_z000.fits \
\
--lmin 2 \
--nbins 30 \
--out_dir ./data/test_box_4/cls \
