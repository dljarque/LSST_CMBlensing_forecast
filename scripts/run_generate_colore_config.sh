#!/bin/bash

# Define the test box to be run
test_box=ensemble_1
/home/dljarque/miniforge3/envs/IFAE/bin/python \
./scripts/generate_colore_config.py \
\
--output ./data/$test_box/ \
--filename param.cfg \
--prefix_out /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/$test_box/results/out \
\
--seed 100000 \
--lensing_nside 1024 \
--z_min 0.0 \
--z_max 3.8 \
\
--omega_M 0.3147 \
--omega_L 0.6853 \
--omega_B 0.04904 \
--h 0.6731 \
--ns 0.9655 \
--As 2.1e-9 \
--w -1.0 \
\
--nz_filenames \
 /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/colore_setup/inputs/nz/nz_lsst_gal.txt \
 /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/colore_setup/inputs/nz/nz_desi_bgs.txt \
 /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/colore_setup/inputs/nz/nz_desi_lrg.txt \
 /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/colore_setup/inputs/nz/nz_desi_elg.txt \
 /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/colore_setup/inputs/nz/nz_desi_qso.txt \
--bz_input_filenames \
 /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/colore_setup/inputs/bias-input/bz-input_lsst_gal.txt \
 /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/colore_setup/inputs/bias-input/bz-input_desi_bgs.txt \
 /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/colore_setup/inputs/bias-input/bz-input_desi_lrg.txt \
 /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/colore_setup/inputs/bias-input/bz-input_desi_elg.txt \
 /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/colore_setup/inputs/bias-input/bz-input_desi_qso.txt \
--Tz_filenames \
 /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/colore_setup/inputs/thresholds/tz_lsst_gal.txt \
 /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/colore_setup/inputs/thresholds/tz_desi_bgs.txt \
 /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/colore_setup/inputs/thresholds/tz_desi_lrg.txt \
 /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/colore_setup/inputs/thresholds/tz_desi_elg.txt \
 /pscratch/sd/d/djarque/LSST_CMBlensing_forecast/data/colore_setup/inputs/thresholds/tz_desi_qso.txt