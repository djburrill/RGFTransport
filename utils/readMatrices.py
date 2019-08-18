# readMatrices.py

# Imports
import numpy as np
import argparse as arg

# Functions
def formatMatrix_QChem(matLines,numBasisFuncs):
    '''
    Format output lines into matrices.
    '''
    # Variables
    matrixRows = []
    rowCounter = 0

    # Intialize matrix
    outMatrix = np.zeros((numBasisFuncs,numBasisFuncs))

    for index,line in enumerate(matLines):
        # Reset counter at headers
        if (index%(numBasisFuncs+1) == 0):
            #print line
            rowCounter = 0
            continue

        # Create new row if does not exist
        if (index < numBasisFuncs+2):
            matrixRows.append([])

        # Read row data
        for matVal in line[1:]:
            matrixRows[rowCounter].append(float(matVal))

        # Increment row counter
        rowCounter += 1

    # Set up matrix
    for indexRow,row in enumerate(matrixRows):
        for indexCol,col in enumerate(row):
            outMatrix[indexRow,indexCol] = col

    return outMatrix

def readMatrices_QChem(outFileName):
    '''
    Read Fock and overlap matrices from QChem output.

    NOTES
        - When running the calculation you should use the rem options 'scf_print 3' and 'IPRINT 200'.
    '''

    # Variables
    FockStartKey = "Extrapolated Alpha Fock Matrix".split()
    FockEndKey = "Alpha MO Eigenvalues".split()
    overlapStartKey = "Overlap Matrix".split()
    overlapEndKey = "Core Hamiltonian".split()
    basisFuncKey = "basis functions".split()
    firstOverlap = False
    numBasisFuncs = 0
    FockRead = False
    overlapRead = False
    FockLineList = []
    overlapLineList = []

    # Read through file
    with open(outFileName,'r') as outFile:
        for line in outFile:
            # Format line
            line = line.strip().split()

            # Skip empty lines
            if (len(line) == 0):
                continue

            # Get number of basis functions
            if (line[-2:] == basisFuncKey):
                #print "READ NUMBER OF BASIS FUNCTIONS"
                numBasisFuncs = int(line[-3])
                continue

            # Check for start keys
            if (line == overlapStartKey):
                #print "START OVERLAP READ"
                if (firstOverlap == False):
                    firstOverlap = True
                else:
                    overlapRead = True
                    overlapLineList = []
                continue

            if (line == FockStartKey):
                #print "START FOCK READ"
                FockRead = True
                FockLineList = []
                continue

            # Check for end keys
            if (line == overlapEndKey):
                #print "STOP OVERLAP READ"
                overlapRead = False
                continue

            if (line == FockEndKey):
                #print "STOP FOCK READ"
                FockRead = False
                continue

            # Read overlap
            if (overlapRead == True):
                overlapLineList.append(line)
                continue

            # Read Fock
            if (FockRead == True):
                FockLineList.append(line)
                continue

    # Format matrices
    overlapMatrix = formatMatrix_QChem(overlapLineList,numBasisFuncs)
    FockMatrix = formatMatrix_QChem(FockLineList,numBasisFuncs)

    # Print matrices to file
    np.savetxt("hamiltonian.dat",FockMatrix)
    np.savetxt("overlap.dat",overlapMatrix)

# Main
if (__name__ == '__main__'):
    # Read command line arguments
    parser = arg.ArgumentParser(description='Grab relevant Hamiltonian and overlap matrices from the output of quantum chemistry calculations.')
    parser.add_argument('-f',
                        '--file',
                        type=str,
                        default='output.out',
                        help='Output file containing Hamiltonian (H), overlap matrices (S).')

    parser.add_argument('-p',
                        '--program',
                        type=str,
                        default='qchem',
                        help='Quantum chemistry program used to calculate properties.')

    args = parser.parse_args()

    # Run program
    if (args.program == 'qchem'):
        readMatrices_QChem(args.file)
