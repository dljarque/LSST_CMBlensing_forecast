import argparse
import sys
from pathlib import Path

# Add the project directory to the Python path
sys.path.append('.')

from utils.colore_param_cfg import (
    build_camb_cosmology,
    save_linear_pk,
    write_colore_param_cfg,
    write_cosmology
)


def main():
    '''
    Generate a CoLoRe parameter config file based on the input cosmological parameters and other inputs.
    '''

    parser = argparse.ArgumentParser()

    # ------------------------------------------------------------------
    # Global arguments
    # ------------------------------------------------------------------
    parser.add_argument("--output", required=True)
    parser.add_argument("--filename", required=True)
    parser.add_argument("--prefix_out", required=True)
    parser.add_argument("--pk_filename", required=False)
    parser.add_argument("--seed", type=int, required=False)
    parser.add_argument("--lensing_nside", type=int, default=1024, required=False)
    parser.add_argument("--z_min", type=float, default=0.0)
    parser.add_argument("--z_max", type=float, default=3.8)


    # ------------------------------------------------------------------
    # Cosmological parameters
    # ------------------------------------------------------------------
    parser.add_argument("--omega_M", type=float, required=True)
    parser.add_argument("--omega_L", type=float, required=True)
    parser.add_argument("--omega_B", type=float, required=True)
    parser.add_argument("--h", type=float, required=True)
    parser.add_argument("--ns", type=float, required=True)
    parser.add_argument("--As", type=float, required=True)
    parser.add_argument("--w", type=float, required=True)

    # ------------------------------------------------------------------
    # Sources
    # ------------------------------------------------------------------
    parser.add_argument("--nz_filenames", nargs='+',  type=str, required=False, default=[], help="List of redshift distribution files.")
    parser.add_argument("--bz_input_filenames", nargs='+', type=str, required=False, default=[], help="List of bias files.")
    parser.add_argument("--Tz_filenames", nargs='+', type=str, required=False, default=[], help="List of threshold files.")

    # ------------------------------------------------------------------
    # Output options
    # ------------------------------------------------------------------
    parser.add_argument("--include_shear", action='store_true')
    parser.add_argument("--include_lensing", action='store_true')
    parser.add_argument("--store_skewers", action='store_true')
    parser.add_argument("--gaussian_skewers", action='store_true')
    parser.add_argument("--output_density", action='store_true')

    # Parse the arguments
    args = parser.parse_args()


    # ------------------------------------------------------------------
    # Prepare output directory
    # ------------------------------------------------------------------

    # Create the output directory if it doesn't exist
    output_dir = Path(args.output)
    prefix_dir = Path(args.prefix_out)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Validate arguments
    # ------------------------------------------------------------------

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

    # Check that the number of nz, bz, and Tz files are the same
    if not (len(args.nz_filenames) == len(args.bz_input_filenames) == len(args.Tz_filenames)):
        parser.error(
            "nz_filenames, bz_input_filenames and Tz_filenames must contain "
            "the same number of files."
        )


    # ------------------------------------------------------------------
    # Build CAMB cosmology
    # ------------------------------------------------------------------


    # Build CAMB cosmology and get the matter power spectrum
    pars, z_star, sigma_8 = build_camb_cosmology(
        omega_M=args.omega_M,
        omega_B=args.omega_B,
        h=args.h,
        As=args.As,
        ns=args.ns,
        w=args.w,
    )


    # ------------------------------------------------------------------
    # Save linear matter power spectrum
    # ------------------------------------------------------------------

    # Get the linear matter power spectrum and save it to a file
    pk_path = save_linear_pk(
        pars,
        output=args.output
    )
    if args.pk_filename is None:
        args.pk_filename = str(prefix_dir.parents[3] / pk_path)


    # ------------------------------------------------------------------
    # Print summary
    # ------------------------------------------------------------------
    print()
    print("=" * 60)
    print("  CoLoRe parameter configuration")
    print("=" * 60)

    print("\nCosmology")
    print(f"  Omega_M : {args.omega_M}")
    print(f"  Omega_L : {args.omega_L}")
    print(f"  Omega_B : {args.omega_B}")
    print(f"  h       : {args.h}")
    print(f"  ns      : {args.ns}")
    print(f"  As      : {args.As:.3e}")
    print(f"  sigma8  : {sigma_8:.5f}")
    print(f"  w       : {args.w}")
    print(f"  z_star  : {z_star:.4f}")

    print("\nSimulation")
    print(f"  seed          : {args.seed}")
    print(f"  lensing nside : {args.lensing_nside}")
    print(f"  z range       : {args.z_min} -> {args.z_max}")

    print("\nOptions")
    print(f"  Shear         : {'yes' if args.include_shear else 'no'}")
    print(f"  Lensing       : {'yes' if args.include_lensing else 'no'}")
    print(f"  Density       : {'yes' if args.output_density else 'no'}")
    print(f"  Skewers       : {'yes' if args.store_skewers else 'no'}")
    print(f"  Gaussian      : {'yes' if args.gaussian_skewers else 'no'}")

    print("\nFiles")
    print(f"  PK            : {str(Path(pk_path))}")
    print(f"  Config        : {output_dir / args.filename}")
    print(f"  Cosmology JSON: {output_dir / 'cosmo_params.json'}")


    # Write the CoLoRe param.cfg file
    write_colore_param_cfg(
        output_file=str(output_dir / args.filename),
        prefix_out=args.prefix_out,
        pk_filename=args.pk_filename,
        nz_filenames=args.nz_filenames,
        bz_filenames=args.bz_input_filenames,
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
        output_file=str(Path(args.output) / "cosmo_params.json"),
        omega_M=args.omega_M,
        omega_L=args.omega_L,
        omega_B=args.omega_B,
        h=args.h,
        As=args.As,
        ns=args.ns,
        w=args.w,
        sigma_8=sigma_8,
        mnu=sum(pars.nu_mass_fractions),
        omk=pars.omk,
        z_out=z_star,
    )

    print("\nConfiguration generated successfully.")
    print("=" * 60)


# Run the main function when the script is executed
if __name__ == "__main__":
    main()