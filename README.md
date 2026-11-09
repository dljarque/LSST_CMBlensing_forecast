# Cross-correlations of CMB lensing with the Lyman alpha forest

At this point just random notes and notebooks, eventually a research project


## How to generate a simulated box

Ask Laura to document the steps here, and upload to GitHub the config file

The simulation output can be found at NERSC:

```/global/cfs/cdirs/desicollab/users/lauracdp/cmb_lensing_box/test_box-0/results```

To see how to read the outputs and make some plots (like CMB lensing x QSOs), checkout: `notebooks/read_colore_outputs.ipynb`


## How to generate the LyaCoLoRe skewers

Ask Laura to document the steps here, and upload to GitHub the config file

The Lya skewers (transmission files) can be found at NERSC:

```/global/cfs/cdirs/desicollab/users/lauracdp/cmb_lensing_box/skewers-0/```


## How to generate the Lya raw deltas

Ask Laura to document the steps here, and upload to GitHub the config file

The raw Lya deltas can be found at NERSC:

```/global/cfs/cdirs/desicollab/users/lauracdp/cmb_lensing_box/mock_analysis/raw_master/deltas_lya/Delta```


## How to rebin the Lya deltas

To see how to read the Lya raw deltas and rebin them into (Nz=20) z bins, checkout: `notebooks/read_raw_deltas.ipynb`

To do this for the entire simulation, one needs to run a script from an interactive node at NERSC: `scripts/rebin_raw_deltas.py`

The rebinned deltas can be found at NERSC under:

```/global/cfs/cdirs/desicollab/users/font/cmb_lensing_fun/rebinned_deltas/```



## How to read and play with the rebinned deltas

Take a look at this notebook to play with them and make Cells plots (also x CMB lensing): `notebooks/plot_lya_2d.ipynb`


This repository started from the history of `cmb_lensing_lya` (original repo by Andreu).
