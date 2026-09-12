#!/bin/bash

/home/dljarque/miniforge3/envs/IFAE/bin/python \
./scripts/generate_rebinned_cls.py 'sources_auto' \
--sources_prefix ./data/test_box_3/results/out_srcs_s1 \
--nside 4096 \
--lmin 2 \
--lmax 2500 \
--zmin_bin 0.0 \
--zmax_bin 3.8 \
--nbins 30 \
--out_dir ./data/test_box_3/cls \
--filename 'desi_qso'