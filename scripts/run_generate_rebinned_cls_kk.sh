#!/bin/bash

/home/dljarque/miniforge3/envs/IFAE/bin/python \
./scripts/generate_rebinned_cls.py 'cmb_lensing_auto' \
--kappa_map ./data/test_box_3/results/out_kappa_z000.fits \
--lmin 2 \
--nbins 30 \
--out_dir ./data/test_box_3/cls \
