# Power spectrum Cls

Here you can find the power spectrum Cls measured from the box.

* The galaxy Cls (`cl_gg`and `cl_kg`) were computed from the sources output of the box. Those were then projected onto a healpy map with resolution:
    * Full gg samples: $\text{NSIDE}=4096$.
    * Full kg samples: $\text{NSIDE}=512$ (matching CoLoRe CMB lensing output map resolution)
    * Tomographic bins sample: $\text{NSIDE}=512$

* The CMB lensing Cls (`cl_kk`) were computed straight from the output map from CoLoRe with resolution $\text{NSIDE}=512$.

* The tomographic redshift bins of some of the power spectra measured are the following:
    * $ \text{Bin0:} \quad  0 < z \leq 2.0 $
    * $ \text{Bin1:} \quad  2.0 < z \leq 3.8 $

