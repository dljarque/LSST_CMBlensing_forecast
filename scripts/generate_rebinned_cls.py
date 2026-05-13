import argparse
from email import parser
import sys
# add the parent directory to the path so we can import from py
sys.path.append('.')

# import the function to save rebinned cls
from py.rebinning import save_rebinned_cls

# main function to run the script
def main():
    '''
    Parse arguments and save rebinned cls for the input box.
    '''
    
    parser = argparse.ArgumentParser()

    # add arguments for the input box and number of bins
    parser.add_argument("--kappa_file", required=True, help="Input CMB lensing convergence (kappa) FITS file")
    parser.add_argument("--out_dir", required=True, help="Output directory for rebinned C_ell products")
    parser.add_argument("--nbins", type=int, default=30, help="Number of bins for rebinning the cls")

    # parse the arguments
    args = parser.parse_args()
    kappa_file = args.kappa_file
    nbins = args.nbins
    out_dir = args.out_dir

    # save the rebinned cls for the input box
    save_rebinned_cls(kappa_file, nbins, out_dir)

    return


if __name__ == "__main__":
    main()

