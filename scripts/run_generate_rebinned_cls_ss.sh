#!/bin/bash

/home/dljarque/miniforge3/envs/IFAE/bin/python \
/home/dljarque/Documents/cmb_lensing_lya/scripts/generate_rebinned_cls.py 'sources_auto' \
--sources 'quasars' \
--sources_folder ./data/third_test_box/quasars/ \
--nside 4096 \
--lmin 2 \
--lmax 2500 \
--zmin_bin 1.5 \
--zmax_bin 3.8 \
--out_dir ./data/third_test_box/quasars/bin1 \
--nbins 30