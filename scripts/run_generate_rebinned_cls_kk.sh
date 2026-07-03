#!/bin/bash

/home/dljarque/miniforge3/envs/IFAE/bin/python \
./scripts/generate_rebinned_cls.py 'cmb_lensing_auto' \
--kappa_map ./data/fourth_test_box/cmb_lensing/out_kappa_z000.fits \
--nbins 30 \
--out_dir ./data/fourth_test_box/cls \
