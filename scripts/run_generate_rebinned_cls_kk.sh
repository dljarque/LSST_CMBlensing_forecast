#!/bin/bash

/home/dljarque/miniforge3/envs/IFAE/bin/python \
/home/dljarque/Documents/cmb_lensing_lya/scripts/generate_rebinned_cls.py 'cmb_lensing_auto' \
--kappa_file ./data/fourth_test_box/cmb_lensing/out_kappa_z000.fits \
--out_dir ./data/fourth_test_box/cmb_lensing/ \
--nbins 30 \
