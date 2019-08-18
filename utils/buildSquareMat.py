# buildSquareMat.py

# Imports
import numpy as np
import matplotlib.pyplot as plt

# Functions
def phi(J,mu,j):
    '''
    Phi function.
    '''

    return np.sqrt(2.0/(J+1))*np.sin(mu*np.pi*j/(J+1.0))

def energy(J,mu,t):
    '''
    Calculate energy weight.
    '''

    return -2*t*np.cos(np.pi*mu/(J+1.0))

def buildSquareMat(t,J):
    '''
    Build Hamiltonian and coupling matrices for square lattice model.
    '''

    # Variables
    hMat = np.zeros([J,J])
    uMat = np.zeros([J,J])

    # Loop over elements in matrices
    # j
    for index1 in range(1,J+1):
        # j'
        for index2 in range(1,J+1):
            # Hamiltonian
            # Sum over mu=1->J
            for mu in range(1,J+1):
                '''
                print "j: " + str(index1)
                print "j': " + str(index2)
                print "mu: " + str(mu)
                print phi(J,mu,index1)
                print phi(J,mu,index2)
                print energy(J,mu,t)
                '''
                hMat[index1-1,index2-1] += phi(J,mu,index1)*phi(J,mu,index2)*energy(J,mu,t)

            # Coupling matrix
            if (index1 == index2):
                uMat[index1-1,index2-1] = -t

    return hMat,uMat

def alpha(t,en,hMat):
    '''
    Build alpha matrix.
    '''

    # Variables
    J = hMat.shape[0]
    alphaOut = np.zeros([J,J])
    enMat = en*np.eye([J,J])

    # Build alpha matrix
    for index1 in range(1,J+1):
        for index2 in range(1,J+1):
            alphaOut[index1,index2] = t*t*(1/(enMat[index1,index2]-hMat[index1,index2]))

    return alphaOut

def epS(hMat,alphaMat):
    '''
    Calculate epsilon^s
    '''

    epS = hMat + alphaMat

    return epS

# Main
if (__name__ == '__main__'):
    # Build square lattice matrix
    t=1
    J=8
    hMat,uMat = buildSquareMat(t,J)

    # Save matrices
    np.savetxt('hamiltonian.dat',hMat)
    np.savetxt('overlap.dat',uMat)
