#!/bin/bash

/home/dljarque/miniforge3/envs/IFAE/bin/python \
/home/dljarque/Documents/cmb_lensing_lya/scripts/generate_rebinned_cls.py 'sources_auto' \
--sources 'quasars' \
--sources_folder ./data/third_test_box/quasars/ \
--nside 4096 \
--lmin 2 \
--lmax 2500 \
--out_dir ./data/third_test_box/quasars/ \
--nbins 30