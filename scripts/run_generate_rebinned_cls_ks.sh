#!/bin/bash

/home/dljarque/miniforge3/envs/IFAE/bin/python \
/home/dljarque/Documents/cmb_lensing_lya/scripts/generate_rebinned_cls.py 'cmblensing_x_sources' \
--kappa_file ./data/third_test_box/cmb_lensing/out_kappa_z000.fits \
--sources 'quasars' \
--sources_folder ./data/third_test_box/quasars/ \
--nside 4096 \
--lmin 2 \
--lmax 2500 \
--zmin_bin 0 \
--zmax_bin 1.5 \
--out_dir ./data/third_test_box/quasars/bin0 \
--nbins 30