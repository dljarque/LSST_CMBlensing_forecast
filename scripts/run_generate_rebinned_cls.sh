#!/bin/bash

/home/dljarque/miniforge3/envs/IFAE/bin/python \
'/home/dljarque/Documents/cmb_lensing_lya/scripts/generate_rebinned_cls.py' \
--kappa_file ./data/third_test_box/out_kappa_z000.fits \
--out_dir ./data/third_test_box \
--nbins 30 \
