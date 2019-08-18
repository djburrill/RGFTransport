# plotTransmission.py

# Imports
import numpy as np
import matplotlib.pyplot as plt

# Functions
def plotTransmission(inFileName):
    '''
    Plot transmission calculated from RGFTranport
    '''
    # Variables
    en = np.arange(-3.0,3.0,0.1)
    realVals = []
    cmplxVals = []
    dataVals = []

    # Load file
    with open(inFileName,'r') as inFile:
        for line in inFile:
            line = line.split(',')
            dataValue = float(line[0][1:])+float(line[1][:-2])*1j
            dataVals.append(np.abs(dataValue))
            realVals.append(float(line[0][1:]))
            cmplxVals.append(float(line[1][:-2]))

    # Plot
    #plt.plot(en,realVals)
    #plt.plot(en,cmplxVals)
    plt.plot(en,dataVals)
    plt.show()

# Main
if (__name__ == '__main__'):
    plotTransmission('output.out')
