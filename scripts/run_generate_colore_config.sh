#!/bin/bash

/home/dljarque/miniforge3/envs/IFAE/bin/python \
'/home/dljarque/Documents/cmb_lensing_lya/scripts/generate_colore_config.py' \
--output ./data/third_test_box/ \
--filename param.cfg \
--prefix_out '$PSCRATCH/cmb_lensing_lya/data/third_test_box/results' \
--omega_M 0.3147 \
--omega_L 0.6853 \
--omega_B 0.04904 \
--h 0.6731 \
--ns 0.9655 \
--As 2.1e-9 \
--w -1.0 \