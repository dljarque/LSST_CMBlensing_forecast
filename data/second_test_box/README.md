# CoLoRe simulation CMB $\kappa$ values

CMB lensing convergence $\kappa$ values of the second test box from the [CoLoRe simulation](https://github.com/damonge/CoLoRe). 
The source plane is the fiducial last scattering surface of the CMB.


## CoLoRe simulation config file

Here you can find the input parameters config file to reproduce this first simulation test box and get the CMB lensing kappa values file.


```python
global =
{
    prefix_out = "/global/cfs/cdirs/desicollab/users/lauracdp/cmb_lensing_box/test_box-no-LSST-0/results/out";
    output_format = "FITS";
    output_density = True;
    pk_filename = "/global/cfs/cdirs/desicollab/users/cramirez/LyaCoLoRe_mocks_inputs/CoLoRe/PlanckDR12_kmax_matterpower_z0.dat";
    z_min = 1.6;
    z_max = 3.8;
    write_pred = False;
    just_write_pred = False;
    seed = 9;
    pred_dz = 0.1;
};

kappa =
{
    z_out = z_max;
    nside = 512;
};

lensing =
{
    nside = 512;
    n_lensing = 11;
    spacing_type = "r";
};

field_par =
{
    r_smooth = 2;
    smooth_potential = True;
    n_grid = 4096;
    dens_type = 0;
    lpt_buffer_fraction = 0.6;
    lpt_interp_type = 1;
    output_lpt = 0;
};

cosmo_par =
{
    omega_M = 0.3147;
    omega_L = 0.6853;
    omega_B = 0.04904;
    h = 0.6731;
    w = -1.0;
    ns = 0.9655;
    sigma_8 = 0.83;
};

srcs1 =
{
    nz_filename = "/global/cfs/cdirs/desicollab/users/lauracdp/photo-z_box/input_files/Nz_qso_130618_2_colore1_hZs.txt";
    bias_filename = "/global/cfs/cdirs/desicollab/users/lauracdp/photo-z_box/input_files/Bz_qso_G18.txt";
    include_shear = False;
    include_lensing = False;
    store_skewers = True;
    gaussian_skewers = True;
};
```