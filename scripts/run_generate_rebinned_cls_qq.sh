#!/bin/bash

/home/dljarque/miniforge3/envs/IFAE/bin/python \
/home/dljarque/Documents/cmb_lensing_lya/scripts/generate_rebinned_cls.py 'quasar_auto' \
--quasars_folder ./data/third_test_box/sources/ \
--out_dir ./data/third_test_box/sources/ \
--nbins 30 \
