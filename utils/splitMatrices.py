# splitMatrices.py

# Imports
import argparse as arg
import numpy as np

# Functions
def splitMatrices(numLeftBasis,numRightBasis,deleteVal):
    '''
    Split matrices to be read by RGFTransport.
    '''
    # Variables

    # Parse delete string
    dList = deleteVal.split(',')

    if (dList[0] == ''):
        dList = []

    for index,dVal in enumerate(dList):
        dList[index] = int(dVal)

    # Read in matrices
    hamiltonian = np.loadtxt("hamiltonian.dat")
    overlap = np.loadtxt("overlap.dat")

    # Remove unwanted rows and columns
    numRemoved = 0

    for dVal in dList:
        hamiltonian = np.delete(hamiltonian,dVal-numRemoved,0)
        hamiltonian = np.delete(hamiltonian,dVal-numRemoved,1)
        overlap = np.delete(overlap,dVal-numRemoved,0)
        overlap = np.delete(overlap,dVal-numRemoved,1)
        numRemoved += 1

    '''
    print "HAMILTONIAN"
    print hamiltonian
    print "\nOVERLAP"
    print overlap
    '''

    # Calculate number of basis functions in device
    numDeviceBasis = hamiltonian.shape[0]-numLeftBasis-numRightBasis

    # Split matrices
    lh = hamiltonian[0:numLeftBasis,0:numLeftBasis]
    ls = overlap[0:numLeftBasis,0:numLeftBasis]
    dh = hamiltonian[numLeftBasis:numLeftBasis+numDeviceBasis,numLeftBasis:numLeftBasis+numDeviceBasis]
    ds = overlap[numLeftBasis:numLeftBasis+numDeviceBasis,numLeftBasis:numLeftBasis+numDeviceBasis]
    rh = hamiltonian[-numRightBasis:,-numRightBasis:]
    rs = overlap[-numRightBasis:,-numRightBasis:]

    # Write matrices
    np.savetxt("lh.dat",lh)
    np.savetxt("ls.dat",ls)
    np.savetxt("dh.dat",dh)
    np.savetxt("ds.dat",ds)
    np.savetxt("rh.dat",rh)
    np.savetxt("rs.dat",rs)

# Main
if (__name__ == '__main__'):
    # Read command line arguments
    parser = arg.ArgumentParser(description='Split matrices to be read by RGFTransport.')

    parser.add_argument('-l',
                        '--numLeftBasis',
                        type=int,
                        help='Number of basis functions in left lead.')

    parser.add_argument('-r',
                        '--numRightBasis',
                        type=int,
                        help='Number of basis functions in right lead.')

    parser.add_argument('-d',
                        '--delete',
                        default='',
                        type=str,
                        help='Remove selected row/columns from matrices (Input comma separated values).')

    args = parser.parse_args()

    # Read, split, and write matrices
    splitMatrices(args.numLeftBasis,args.numRightBasis,args.delete)
