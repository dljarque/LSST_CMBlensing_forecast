import argparse
import sys
# add the parent directory to the path so we can import from py
sys.path.append('.')

from py.rebinning import save_rebinned_cls_kk, save_rebinned_cls_qq

def cmb_lensing_auto(args):
    '''Saves rebinned cls for CMB lensing.'''
    save_rebinned_cls_kk(args.kappa_file, args.nbins, args.out_dir)

def quasar_auto(args):
    '''Saves rebinned cls for Quasars.'''
    save_rebinned_cls_qq(args.quasars_folder, args.nbins, args.out_dir)


if __name__ == "__main__":
    # Set up the main parser
    parser = argparse.ArgumentParser(description="Rebinning Cls Pipeline")
    subparsers = parser.add_subparsers(dest="task", required=True, help="Task to run")

    # Create the 'cmb_lensing_auto' sub-command
    parser_cmb = subparsers.add_parser("cmb_lensing_auto", help="Run CMB lensing auto-correlation")
    parser_cmb.add_argument("--kappa_file", required=True, help="Input CMB lensing convergence FITS file")
    parser_cmb.add_argument("--out_dir", required=True, help="Output directory")
    parser_cmb.add_argument("--nbins", type=int, default=30, help="Number of bins")
    parser_cmb.set_defaults(func=cmb_lensing_auto) # Link to the function

    # Create the 'quasar_auto' sub-command
    parser_quasar = subparsers.add_parser("quasar_auto", help="Run Quasar auto-correlation")
    parser_quasar.add_argument("--quasars_folder", required=True, help="Input folder with quasar data")
    parser_quasar.add_argument("--out_dir", required=True, help="Output directory")
    parser_quasar.add_argument("--nbins", type=int, default=30, help="Number of bins")
    parser_quasar.set_defaults(func=quasar_auto) # Link to the function

    # Parse and execute automatically
    args = parser.parse_args()
    args.func(args) # This calls the linked function and passes the args to it