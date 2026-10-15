#!/bin/bash

PYTHON=/home/dljarque/miniforge3/envs/IFAE/bin/python
SCRIPT=./scripts/compute_colore_cls.py

# Set the test box, sources prefix, and output filename
test_box='test_box_5'

# CMB lensing auto-spectrum
$PYTHON $SCRIPT 'cmb_lensing_auto' \
--kappa_map ./data/$test_box/results/out_kappa_z000.fits \
--lmin 2 \
--nbins 30 \
--out_dir ./data/$test_box/cls \
