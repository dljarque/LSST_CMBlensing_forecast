import numpy as np
import camb
import json


def build_camb_cosmology(
    omega_M,
    omega_B,
    h,
    ns,
    As,
    w=-1.0,
):
    """
    Create CAMB cosmology and compute z_star.
    """

    H0 = 100 * h
    ombh2 = omega_B * h**2
    omch2 = (omega_M - omega_B) * h**2

    pars = camb.CAMBparams()

    pars.set_cosmology(
        H0=H0,
        ombh2=ombh2,
        omch2=omch2,
    )

    pars.InitPower.set_params(
        ns=ns,
        As=As,  # temporary initial guess
    )

    pars.set_dark_energy(w=w)

    # results object
    pars.WantTransfer = True
    results = camb.get_results(pars)

    # recombination redshift
    z_star = results.get_derived_params()["zstar"]
    sigma_8 = results.get_sigma8_0()

    return pars, z_star, sigma_8


def save_linear_pk(
    pars,
    output,
    kmax=10.0,
    nk=500,
):
    """
    Compute and save linear P(k, z=0).
    """

    # Set the matter power spectrum parameters
    pars.set_matter_power(
        redshifts=[0.0],
        kmax=kmax
    )

    results = camb.get_results(pars)

    # Get P(k) at z=0
    k, z, pk = results.get_matter_power_spectrum(
        minkh=1e-4,
        maxkh=kmax,
        npoints=nk
    )
    
    # Extract P(k) at z=0 and save to file
    pk0 = pk[0]
    # Combine k and P(k) into a single array for saving
    data = np.column_stack([k, pk0])

    # Save to file
    output_pk_file = f"{output}pk_linear.dat"
    np.savetxt(output_pk_file, data)

    print(f"Saved P(k) to {output_pk_file}")

    return output_pk_file


def write_colore_param_cfg(
    output_file,
    prefix_out,
    pk_filename,
    nz_filename,
    bias_filename,
    z_out,
    omega_M,
    omega_L,
    omega_B,
    h,
    w,
    sigma_8,
    z_min=1.6,
    z_max=3.8,
    seed=9,
    nside=512,

):

    """
    Generate CoLoRe param.cfg file from cosmological parameters and other inputs.

    Parameters
    ----------
    output_file : str
        Path where param.cfg will be written
    prefix_out : str
        Output prefix for CoLoRe simulationproducts
    pk_filename : str
        Input matter power spectrum file
    omega_M, omega_L, omega_B, h, w, sigma_8 : float
        Cosmological parameters
    nz_filename : str, optional
        Redshift distribution file (default: input/nz.txt)
    bias_filename : str, optional
        Galaxy bias file (default: input/bias.txt)
    """

    content = f"""
global =
{{
    prefix_out = "{prefix_out}";
    output_format = "FITS";
    output_density = True;
    pk_filename = "{pk_filename}";
    z_min = {z_min};
    z_max = {z_max};
    write_pred = False;
    just_write_pred = False;
    seed = {seed};
    pred_dz = 0.1;
}};

kappa =
{{
    z_out = {z_out:.2f};
    nside = {nside};
}};

lensing =
{{
    nside = {nside};
    n_lensing = 11;
    spacing_type = "r";
}};

field_par =
{{
    r_smooth = 2;
    smooth_potential = True;
    n_grid = 4096;
    dens_type = 0;
    lpt_buffer_fraction = 0.6;
    lpt_interp_type = 1;
    output_lpt = 0;
}};

cosmo_par =
{{
    omega_M = {omega_M};
    omega_L = {omega_L};
    omega_B = {omega_B};
    h = {h};
    w = {w};
    sigma_8 = {sigma_8:.2f};
}};

srcs1 =
{{
    nz_filename = "{nz_filename}";
    bias_filename = "{bias_filename}";
    include_shear = False;
    include_lensing = False;
    store_skewers = True;
    gaussian_skewers = True;
}};
"""

    with open(output_file, "w") as f:
        f.write(content.strip())

    print(f"Saved config to {output_file}")


def write_cosmology(
    output_file,
    omega_M,
    omega_L,
    omega_B,
    h,
    ns,
    As,
    w,
    sigma_8,
    z_out,
    z_max=3.8,
):

    """
    Generate a JSON file with cosmological parameters and other CoLoRe params for record-keeping and reproducibility.
    """

    cosmo_par = {
        "omega_M": omega_M,
        "omega_L": omega_L,
        "omega_B": omega_B,
        "h": h,
        "As": As,
        "ns": ns,
        "w": w,
        "sigma_8": round(float(sigma_8), 4),
        "z_max": z_max,
        "z_out": round(float(z_out), 2),
    }

    with open(output_file, "w") as f:
        json.dump(cosmo_par, f, indent=4)

    print(f"Saved cosmo params to {output_file}")