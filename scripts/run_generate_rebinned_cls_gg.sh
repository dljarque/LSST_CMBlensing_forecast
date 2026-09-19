#!/bin/bash

/home/dljarque/miniforge3/envs/IFAE/bin/python \
./scripts/generate_rebinned_cls.py 'sources_auto' \
--sources_prefix ./data/test_box_4/results/out_srcs_s2 \
--nside 1024 \
--lmin 2 \
--lmax 2500 \
--nbins 30 \
--out_dir ./data/test_box_4/cls \
--filename 'desi_bgs' \
#--zmean 2.0 \
#--dz 0.2 \
#--sigma_z 0.05 \
