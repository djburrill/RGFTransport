# plotSGF.py

# Imports
import matplotlib.pyplot as plt
import numpy as np

# Functions
def plotSGF(inFileName):
    '''
    Plot SGF to compare to literature
    '''

    # Variables

    # Load file
    energy,real1,imag1 = np.loadtxt(inFileName,unpack=True)

    # Plot data
    plt.plot(energy,real1)
    plt.title("Surface Green's Function (Real)")
    plt.show()

    plt.plot(energy,imag1)
    plt.title("Surface Green's Function (Imaginary)")
    plt.show()

def plotSGF_v1(inFileName):
    '''
    Plot SGF to compare to literature
    '''

    # Variables

    # Load file
    real1,imag1,real2,imag2,real3,imag3,real4,imag4 = np.loadtxt(inFileName,unpack=True)

    # Energies
    en = np.arange(-6.0,6.0,0.001)

    # Plot data
    plt.plot(en,real1)
    plt.plot(en,real2)
    plt.plot(en,real3)
    plt.plot(en,real4)
    plt.show()

    plt.plot(en,imag1)
    plt.plot(en,imag2)
    plt.plot(en,imag3)
    plt.plot(en,imag4)
    plt.show()

# Main
if (__name__ == '__main__'):
    plotSGF('sgf.dat')
