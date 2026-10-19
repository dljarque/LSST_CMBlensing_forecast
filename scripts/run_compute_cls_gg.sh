#!/bin/bash

PYTHON=/home/dljarque/miniforge3/envs/IFAE/bin/python
SCRIPT=./scripts/compute_colore_cls.py

# Set the test box, sources prefix, and output filename
test_box='test_box_6_old'

# Galaxy auto-spectrum
$PYTHON $SCRIPT 'sources_auto' \
    --sources_prefix ./data/$test_box/results/out_srcs_s5 \
    --nside 1024 \
    --lmin 2 \
    --lmax 1024 \
    --nbins 30 \
    --out_dir ./data/$test_box/cls \
    --filename 'desi_qso_bin0' \
    --zmean 2.0 \
    --dz 0.2 \
    --sigma_z 0.05
