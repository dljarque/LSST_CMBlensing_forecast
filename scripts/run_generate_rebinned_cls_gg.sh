#!/bin/bash

/home/dljarque/miniforge3/envs/IFAE/bin/python \
./scripts/generate_rebinned_cls.py 'sources_auto' \
--sources_prefix ./data/fourth_test_box/results/out_srcs_s5 \
--nside 512 \
--lmin 2 \
--zmin_bin 2.0 \
--zmax_bin 3.8 \
--nbins 30 \
--out_dir ./data/fourth_test_box/cls \
--filename 'desi_qso_bin1'