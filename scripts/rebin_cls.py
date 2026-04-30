import argparse
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
    parser.add_argument("--input", required=True, help="Input data box")
    parser.add_argument("--nbins", type=int, default=40, help="Number of bins for rebinning the cls")

    # parse the arguments
    args = parser.parse_args()
    data = args.input
    nbins = args.nbins

    # save the rebinned cls for the input box
    save_rebinned_cls(data, nbins)

    return


if __name__ == "__main__":
    main()

