import argparse
import sys
# add the parent directory to the path so we can import from py
sys.path.append('.')

from utils.colore_param_cfg import *


def main():
    '''
    Generate a CoLoRe parameter config file based on the input cosmological parameters and other inputs.
    '''

    parser = argparse.ArgumentParser()

    # Global arguments
    parser.add_argument("--output", required=True)
    parser.add_argument("--filename", required=True)
    parser.add_argument("--prefix_out", required=True)
    parser.add_argument("--pk_filename", required=False)
    parser.add_argument("--seed", type=int, required=False)
    parser.add_argument("--lensing_nside", type=int, required=False)
    parser.add_argument("--z_min", type=float, default=1.6)
    parser.add_argument("--z_max", type=float, default=3.8)


    # Cosmological parameters
    parser.add_argument("--omega_M", type=float, required=True)
    parser.add_argument("--omega_L", type=float, required=True)
    parser.add_argument("--omega_B", type=float, required=True)
    parser.add_argument("--h", type=float, required=True)
    parser.add_argument("--ns", type=float, required=True)
    parser.add_argument("--As", type=float, default=2e-9)
    parser.add_argument("--w", type=float, default=-1.0)

    # Sources
    parser.add_argument("--nz_filenames", nargs='+',  type=str, required=False, default=[], help="List of redshift distribution files.")
    parser.add_argument("--bz_filenames", nargs='+', type=str, required=False, default=[], help="List of bias files.")
    parser.add_argument("--Tz_filenames", nargs='+', type=str, required=False, default=[], help="List of threshold files.")
    parser.add_argument("--include_shear", action='store_true')
    parser.add_argument("--include_lensing", action='store_true')
    parser.add_argument("--store_skewers", action='store_true')
    parser.add_argument("--gaussian_skewers", action='store_true')
    parser.add_argument("--output_density", action='store_true')

    # Parse the arguments
    args = parser.parse_args()

    # Print the input parameters for verification
    print("Shear is included") if args.include_shear else print("Shear is NOT included")
    print("Lensing is included") if args.include_lensing else print("Lensing is NOT included")
    if args.store_skewers:
        if args.gaussian_skewers:
            print("Skewers will be stored and will be Gaussian")
        else:
            print("Skewers will be stored but will NOT be Gaussian")
    else:
        print("Skewers will NOT be stored")
    print("Density field will be stored") if args.output_density else print("Density field will NOT be stored")


    # Build CAMB cosmology and get the matter power spectrum
    pars, z_star, sigma_8 = build_camb_cosmology(
        omega_M=args.omega_M,
        omega_B=args.omega_B,
        h=args.h,
        As=args.As,
        ns=args.ns,
        w=args.w,
    )


    # Get the linear matter power spectrum and save it to a file
    pk_filename = save_linear_pk(
        pars,
        output=args.output
    )
    if args.pk_filename is None:
        args.pk_filename = f"/pscratch/sd/d/djarque/cmb_lensing_lya/{pk_filename}".replace("/./", "/")  # Ensure no double slashes in the path


    # Write the CoLoRe param.cfg file
    write_colore_param_cfg(
        output_file=args.output+args.filename,
        prefix_out=args.prefix_out,
        pk_filename=args.pk_filename,
        nz_filenames=args.nz_filenames,
        bz_filenames=args.bz_filenames,
        Tz_filenames=args.Tz_filenames,
        seed=args.seed,
        nside=args.lensing_nside,
        z_min=args.z_min,
        z_max=args.z_max,
        z_out=z_star,
        omega_M=args.omega_M,
        omega_L=args.omega_L,
        omega_B=args.omega_B,
        h=args.h,
        w=args.w,
        ns=args.ns,
        sigma_8=sigma_8,
        include_shear=args.include_shear,
        include_lensing=args.include_lensing,
        store_skewers=args.store_skewers,
        gaussian_skewers=args.gaussian_skewers,
        output_density=args.output_density,
    )


    # Write the cosmology.json file for record-keeping
    write_cosmology(
        output_file=args.output+"cosmo_params.json",
        omega_M=args.omega_M,
        omega_L=args.omega_L,
        omega_B=args.omega_B,
        h=args.h,
        As=args.As,
        ns=args.ns,
        w=args.w,
        sigma_8=sigma_8,
        z_out=z_star,
    )


# Run the main function when the script is executed
if __name__ == "__main__":
    main()