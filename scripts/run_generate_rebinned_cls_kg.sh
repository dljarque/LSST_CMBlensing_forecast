#!/bin/bash

/home/dljarque/miniforge3/envs/IFAE/bin/python \
./scripts/generate_rebinned_cls.py 'cmblensing_x_sources' \
--kappa_map ./data/test_box_3/results/out_kappa_z000.fits \
--sources_prefix ./data/test_box_3/results/out_srcs_s1 \
--lmin 2 \
--zmin_bin 0.0 \
--zmax_bin 3.8 \
--nbins 30 \
--out_dir ./data/test_box_3/cls \
--filename 'desi_qso'