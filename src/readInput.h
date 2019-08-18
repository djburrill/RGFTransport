/** Read input files.
*/
#ifndef _readInput_h_
#define _readInput_h_
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <atom.h>
#include <armadillo>

/** Class to read input files
* Input files will have a specific format. This class parses the format and saves all of the parameter information.
*/
class readInput {
public:
  // Public variables
  double m_minSliceDist;        /**< Minimum distance between slices */
  int m_llucAtoms;              /**< Number of atoms per unit cell in left lead */
  int m_rlucAtoms;              /**< Number of atoms per unit cell in right lead */
  int m_numBasisUCL = 0;        /**< Numer of basis functions per unit cell in left lead */
  int m_numBasisUCR = 0;        /**< Numer of basis functions per unit cell in right lead */
  double m_imag_eta = 0.000001; /**< Imaginary component of energy */
  double m_sc_conv = 0.000001;  /**< Self-consistent cycle convergence tolderance */
  double m_vlow = -1.0;         /**< Low voltage */
  double m_vhigh = 1.0;         /**< High voltage */
  double m_dv = 0.1;            /**< Voltage spacing */
  double m_elow = -2.0;         /**< Low energy (from energy integral) */
  double m_ehigh = 2.0;         /**< High energy (from energy integral) */
  double m_de = 0.001;          /**< Energy spacing (from energy integral) */

  /** Constructor
  */
  readInput(std::string a_inFileName);

  /** Get device molecule
  */
  std::vector<atom> getMolecule(){return m_molecule;};

  /** Get left lead atoms
  */
  std::vector<atom> getleftLead(){return m_lLead;};

  /** Get right lead atoms
  */
  std::vector<atom> getRightLead(){return m_rLead;};

  /** Trim string from left of specified character
  */
  std::string trimStringLeft(std::string inString, char trimChar);

  /** Trim string from left of spaces
  */
  std::string trimStringLeft(std::string inString);

  /** Trim string from right of specified character
  */
  std::string trimStringRight(std::string inString, char trimChar);

  /** Trim string from right of spaces
  */
  std::string trimStringRight(std::string inString);

  /** Trim string on both sides of specified character
  */
  std::string trimString(std::string inString, char trimChar);

  /** Trim string on both sides of spaces
  */
  std::string trimString(std::string inString);

  /** Print basis numbers
  */
  void print_Basis();

  /** Read matrix from file
  */
  void readMatrix(arma::mat& a_inMat, std::string a_fileName);

  /** Get left lead Hamiltonian matrix
  */
  const arma::mat& getLHMat(){return m_lhMat;};

  /** Get left lead overlap matrix
  */
  const arma::mat& getLSMat(){return m_lsMat;};

  /** Get right lead Hamiltonian matrix
  */
  const arma::mat& getRHMat(){return m_rhMat;};

  /** Get right lead overlap matrix
  */
  const arma::mat& getRSMat(){return m_rsMat;};

  /** Get device Hamiltonian matrix
  */
  const arma::mat& getDHMat(){return m_dhMat;};

  /** Get device overlap matrix
  */
  const arma::mat& getDSMat(){return m_dsMat;};

  /** Calculate number of basis functions per unit cell in the leads
  */
  void calcBasisUCLead();

private:
  // Private variables
  std::string m_inputCard = "None";
  std::vector<atom> m_molecule;         /**< Vector of atoms in device */
  std::vector<atom> m_lLead;            /**< Vector of atoms in left lead */
  std::vector<atom> m_rLead;            /**< Vector of atoms in right lead */
  std::map<std::string,int> m_numBasis; /**< Number of basis functions per atomic label */

  // Matrices
  arma::mat m_lhMat;    /**< Left lead Hamiltonian matrix */
  arma::mat m_lsMat;    /**< Left lead overlap matrix */
  arma::mat m_rhMat;    /**< Right lead Hamiltonian matrix */
  arma::mat m_rsMat;    /**< Right lead overlap matrix */
  arma::mat m_dhMat;    /**< Device Hamiltonian matrix */
  arma::mat m_dsMat;    /**< Device overlap matrix */

  /** Read input file
  */
  void readInputFile(std::string a_inFileName);

  /** Read property variables
  */
  void readProperty(std::string a_inString);

  /** Read atomic positions
  */
  void readAtomPos(std::vector<atom>& a_molecule, std::string a_inString);

  /** Read number of basis functions
  */
  void readBasis(std::string a_inString);

};
#endif
