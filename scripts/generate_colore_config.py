import argparse
import sys
# add the parent directory to the path so we can import from py
sys.path.append('.')

from py.colore_param_cfg import *


def main():
    '''
    Generate a CoLoRe parameter config file based on the input cosmological parameters and other inputs.
    '''

    parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument("--output", required=True)
    parser.add_argument("--filename", required=True)
    parser.add_argument("--prefix_out", required=True)
    parser.add_argument("--pk_filename", required=False)
    parser.add_argument("--nz_filename",
                        default="pscratch/sd/d/djarque/cmb_lensing_lya/data/Nz_qso_130618_2_colore1_hZs.txt")
    parser.add_argument("--bias_filename",
                        default="pscratch/sd/d/djarque/cmb_lensing_lya/data/Bz_qso_G18.txt")

    # Cosmological parameters
    parser.add_argument("--omega_M", type=float, required=True)
    parser.add_argument("--omega_L", type=float, required=True)
    parser.add_argument("--omega_B", type=float, required=True)
    parser.add_argument("--h", type=float, required=True)
    parser.add_argument("--ns", type=float, required=True)
    parser.add_argument("--As", type=float, default=2e-9)
    parser.add_argument("--w", type=float, default=-1.0)

    args = parser.parse_args()


    # Build CAMB cosmology and get the matter power spectrum
    pars, z_star = build_camb_cosmology(
        omega_M=args.omega_M,
        omega_B=args.omega_B,
        h=args.h,
        ns=args.ns,
        As=args.As,
        w=args.w,
    )


    # Get the linear matter power spectrum and save it to a file
    pk_filename = save_linear_pk(
        pars,
        output=args.output
    )
    if args.pk_filename is None:
        args.pk_filename = f"pscratch/sd/d/djarque/cmb_lensing_lya/{pk_filename}"

    # Write the CoLoRe param.cfg file
    write_colore_param_cfg(
        output_file=args.output+args.filename,
        prefix_out=args.prefix_out,
        pk_filename=pk_filename,
        nz_filename=args.nz_filename,
        bias_filename=args.bias_filename,
        z_out = z_star,
        omega_M=args.omega_M,
        omega_L=args.omega_L,
        omega_B=args.omega_B,
        h=args.h,
        ns=args.ns,
        As=args.As,
        w=args.w
    )


    # Write the cosmology.json file for record-keeping
    write_cosmology(
        output_file=args.output+"cosmo_params.json",
        omega_M=args.omega_M,
        omega_L=args.omega_L,
        omega_B=args.omega_B,
        h=args.h,
        ns=args.ns,
        As=args.As,
        w=args.w,
        z_out = z_star,
    )


# Run the main function when the script is executed
if __name__ == "__main__":
    main()