import argparse
import sys
# add the parent directory to the path so we can import from py
sys.path.append('.')

from utils.rebinning import save_rebinned_cls_kk, save_rebinned_cls_ss, save_rebinned_cls_ks

def cmb_lensing_auto(args):
    '''Saves rebinned Cls for CMB lensing.'''
    save_rebinned_cls_kk(args.kappa_map, args.lmin, args.lmax, args.nbins, args.out_dir)

def sources_auto(args):
    '''Saves rebinned Cls for Sources.'''
    save_rebinned_cls_ss(args.sources_prefix, args.nside, args.lmin, args.lmax, args.zmin_bin, args.zmax_bin, args.nbins, args.out_dir, args.filename)

def cmblensing_x_sources(args):
    '''Saves rebinned Cls for CMB lensing x Sources.'''
    save_rebinned_cls_ks(args.kappa_map, args.sources_prefix, args.nside, args.lmin, args.lmax, args.zmin_bin, args.zmax_bin, args.nbins, args.out_dir, args.filename)


if __name__ == "__main__":
    # Set up the main parser
    parser = argparse.ArgumentParser(description="Rebinning Cls Pipeline")
    subparsers = parser.add_subparsers(dest="task", required=True, help="Task to run")

    # Create the 'cmb_lensing_auto' sub-command
    parser_cmb = subparsers.add_parser("cmb_lensing_auto", help="Run CMB lensing auto-correlation")
    parser_cmb.add_argument("--kappa_map", required=True, help="Input CMB lensing convergence FITS file")
    parser_cmb.add_argument("--lmin", type=int, default=None, help="Minimum multipole moment to include")
    parser_cmb.add_argument("--lmax", type=int, default=None, help="Maximum multipole moment to include")
    parser_cmb.add_argument("--nbins", type=int, default=30, help="Number of bins")
    parser_cmb.add_argument("--out_dir", required=True, help="Output directory")
    parser_cmb.set_defaults(func=cmb_lensing_auto) # Link to the function

    # Create the 'sources_auto' sub-command
    parser_source = subparsers.add_parser("sources_auto", help="Run Sources auto-correlation")
    parser_source.add_argument("--sources_prefix", required=True, help="Input folder with sources data")
    parser_source.add_argument("--nside", type=int, required=True, help="Healpix resolution for the delta_g map (must be a power of 2)")
    parser_source.add_argument("--lmin", type=int, default=None, help="Minimum multipole moment to include")
    parser_source.add_argument("--lmax", type=int, default=None, help="Maximum multipole moment to include")
    parser_source.add_argument("--zmin_bin", type=float, default=None, help="Lower limit of the tomographic bin")
    parser_source.add_argument("--zmax_bin", type=float, default=None, help="Upper limit of the tomographic bin")
    parser_source.add_argument("--nbins", type=int, default=30, help="Number of bins")
    parser_source.add_argument("--out_dir", required=True, help="Output directory")
    parser_source.add_argument("--filename", required=True, help="Output cls filename")
    parser_source.set_defaults(func=sources_auto) # Link to the function

    # Create the 'cmblensing_x_sources' sub-command
    parser_cmbxsource = subparsers.add_parser("cmblensing_x_sources", help="Run CMB lensing x sources cross-correlation")
    parser_cmbxsource.add_argument("--kappa_map", required=True, help="Input CMB lensing convergence FITS file")
    parser_cmbxsource.add_argument("--sources_prefix", required=True, help="Input folder with sources data")
    parser_cmbxsource.add_argument("--nside", type=int, required=True, help="Healpix resolution for the delta_g map (must be a power of 2)")
    parser_cmbxsource.add_argument("--lmin", type=int, default=None, help="Minimum multipole moment to include")
    parser_cmbxsource.add_argument("--lmax", type=int, default=None, help="Maximum multipole moment to include")
    parser_cmbxsource.add_argument("--zmin_bin", type=float, default=None, help="Lower limit of the tomographic bin")
    parser_cmbxsource.add_argument("--zmax_bin", type=float, default=None, help="Upper limit of the tomographic bin")
    parser_cmbxsource.add_argument("--nbins", type=int, default=30, help="Number of bins")
    parser_cmbxsource.add_argument("--out_dir", required=True, help="Output directory")
    parser_cmbxsource.add_argument("--filename", required=True, help="Output cls filename")
    parser_cmbxsource.set_defaults(func=cmblensing_x_sources) # Link to the function

    # Parse and execute automatically
    args = parser.parse_args()
    args.func(args) # This calls the linked function and passes the args to it