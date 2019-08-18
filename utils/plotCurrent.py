# plotCurrent.py

# Imports
import numpy as np
import matplotlib.pyplot as plt

# Functions
def plotDOS(inFileName):
    '''
    Plot density of states
    '''

    # Load file
    enValList,dataVals = np.loadtxt(inFileName,unpack=True)

    # Plot DOS
    plt.plot(enValList,dataVals)
    plt.title('Density of States')
    plt.show()

def plotTransmission(inFileName):
    '''
    Plot transmission calculated from RGFTranport
    '''
    # Variables
    en = np.arange(-3.0,3.0,0.1)
    realVals = []
    cmplxVals = []
    dataVals = []
    enValList = []
    condList = []

    # Load file
    '''
    with open(inFileName,'r') as inFile:
        for line in inFile:
            line = line.split()
            line1 = line[1].split(',')
            dataValue = float(line1[0][1:])+float(line1[1][:-1])*1j
            dataVals.append(np.abs(dataValue))
            realVals.append(float(line1[0][1:]))
            cmplxVals.append(float(line1[1][:-2]))
            enValList.append(float(line[0]))
    '''
    enValList,dataVals = np.loadtxt("IV.dat",unpack=True)

    '''
    # Adjust abs data for voltage
    for index,datVal in enumerate(dataVals):
        if (enValList[index] < 0.0):
            dataVals[index] = -datVal
    '''

    # Plot
    #plt.plot(en,realVals)
    #plt.plot(en,cmplxVals)
    plt.plot(enValList,dataVals)
    plt.title('Current vs. Voltage')
    plt.show()

    # Calculate conductivity
    for index in range(len(dataVals)-1):
        condList.append((dataVals[index+1]-dataVals[index])/(enValList[index+1]-enValList[index]))

    plt.plot(enValList[:-1],condList)
    plt.title('Conductivity vs. Voltage')
    plt.show()

# Main
if (__name__ == '__main__'):
    plotTransmission('output.out')
    plotDOS("DOS.dat")
