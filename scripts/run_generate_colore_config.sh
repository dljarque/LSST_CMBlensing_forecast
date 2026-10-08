#!/bin/bash

/home/dljarque/miniforge3/envs/IFAE/bin/python \
./scripts/generate_colore_config.py \
\
--output ./data/test_box_6/ \
--filename param.cfg \
--prefix_out /pscratch/sd/d/djarque/cmb_lensing_lya/data/test_box_6/results/out \
\
--seed 66 \
--lensing_nside 1024 \
--z_min 0 \
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
 /pscratch/sd/d/djarque/cmb_lensing_lya/data/colore_files/nz/nz_lsst_gal.txt \
 /pscratch/sd/d/djarque/cmb_lensing_lya/data/colore_files/nz/nz_desi_bgs.txt \
 /pscratch/sd/d/djarque/cmb_lensing_lya/data/colore_files/nz/nz_desi_lrg.txt \
 /pscratch/sd/d/djarque/cmb_lensing_lya/data/colore_files/nz/nz_desi_elg.txt \
 /pscratch/sd/d/djarque/cmb_lensing_lya/data/colore_files/nz/nz_desi_qso.txt \
--bz_filenames \
 /pscratch/sd/d/djarque/cmb_lensing_lya/data/colore_files/bz-tunning/bz-tunning_lsst_gal.txt \
 /pscratch/sd/d/djarque/cmb_lensing_lya/data/colore_files/bz-tunning/bz-tunning_desi_bgs.txt \
 /pscratch/sd/d/djarque/cmb_lensing_lya/data/colore_files/bz-tunning/bz-tunning_desi_lrg.txt \
 /pscratch/sd/d/djarque/cmb_lensing_lya/data/colore_files/bz-tunning/bz-tunning_desi_elg.txt \
 /pscratch/sd/d/djarque/cmb_lensing_lya/data/colore_files/bz-tunning/bz-tunning_desi_qso.txt \
--Tz_filenames \
 /pscratch/sd/d/djarque/cmb_lensing_lya/data/colore_files/Tz/Tz_lsst_gal.txt \
 /pscratch/sd/d/djarque/cmb_lensing_lya/data/colore_files/Tz/Tz_desi_bgs.txt \
 /pscratch/sd/d/djarque/cmb_lensing_lya/data/colore_files/Tz/Tz_desi_lrg.txt \
 /pscratch/sd/d/djarque/cmb_lensing_lya/data/colore_files/Tz/Tz_desi_elg.txt \
 /pscratch/sd/d/djarque/cmb_lensing_lya/data/colore_files/Tz/Tz_desi_qso.txt