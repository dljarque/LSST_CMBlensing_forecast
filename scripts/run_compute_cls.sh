#!/bin/bash

PYTHON=/home/dljarque/miniforge3/envs/IFAE/bin/python
SCRIPT=./scripts/compute_colore_cls.py

# Set the test box, sources prefix, and output filename
test_box='ensemble_10'
sources='out_srcs_s5'
filename='desi_qso_bin0'

# Set to true to run tomographic binning and include zmean, dz, and sigma_z. Set to false to exclude them.
TOMOGRAPHIC_BINNING=true

zmean=2.0
dz=0.2
sigma_z=0.05

Z_ARGS=""
if [ "$TOMOGRAPHIC_BINNING" = true ]; then
    Z_ARGS="--zmean $zmean --dz $dz --sigma_z $sigma_z"
fi


# Galaxy auto-spectrum
$PYTHON $SCRIPT 'sources_auto' \
    --sources_prefix ./data/$test_box/results/$sources \
    --nside 1024 \
    --lmin 2 \
    --nbins 30 \
    --out_dir ./data/$test_box/cls \
    --filename $filename \
    $Z_ARGS

# CMB lensing x galaxy cross-spectrum
$PYTHON $SCRIPT 'cmblensing_x_sources' \
    --kappa_map ./data/$test_box/results/out_kappa_z000.fits \
    --sources_prefix ./data/$test_box/results/$sources \
    --lmin 2 \
    --nbins 30 \
    --out_dir ./data/$test_box/cls \
    --filename $filename \
    $Z_ARGS