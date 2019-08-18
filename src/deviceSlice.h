/** Partition device into slices.
*/
#ifndef _deviceSlice_h_
#define _deviceSlice_h_
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <armadillo>
#include <atom.h>

class deviceSlice {
public:
  // Public functions

  /** Constructor
  */
  deviceSlice(const std::vector<atom>& a_molecule,
              const arma::mat& a_HMat,
              const arma::mat& a_SMat,
              double a_minSliceDist,
              int a_llucAtoms,
              int a_rlucAtoms);

  /** Slice molecule
    * Transport will be along the x direction so slices will be taken along this axis.
  */
  void slice();

  /** Slice Hamiltonian and overlap matrices
    * Slice up the device Hamiltonian and overlap matrices into smaller Hamiltonians and the coupling matrices. uSlices should be one element shorter than HamSlices.
  */
  void sliceHam(std::vector<arma::cx_mat>& a_HamSlices,
                std::vector<arma::cx_mat>& a_overSlices,
                std::vector<arma::cx_mat>& a_uSlices,
                std::vector<arma::cx_mat>& a_suSlices,
                std::vector<std::vector<atom> >& a_slices,
                arma::mat& a_HMat,
                arma::mat& a_SMat);

  /** Get Hamiltonian slices
  */
  const std::vector<arma::cx_mat> getHamSlices(){return m_HamSlices;};

  /** Get overlap matrix slices
  */
  const std::vector<arma::cx_mat> getOverSlices(){return m_overSlices;};

  /** Get Hamiltonian coupling slices
  */
  const std::vector<arma::cx_mat> getUSlices(){return m_uSlices;};

  /** Get overlap coupling slices
  */
  const std::vector<arma::cx_mat> getSUSlices(){return m_suSlices;};

private:
  std::vector<atom> m_molecule;                     /**< Vector of atoms in molecule*/
  std::vector<std::vector<atom> > m_slices;         /**< Vector of vector of atoms in slices*/
  std::vector<arma::cx_mat> m_HamSlices;            /**< Vector of Hamiltonians of slices*/
  std::vector<arma::cx_mat> m_overSlices;           /**< Vector of overlap matrices of slices*/
  std::vector<arma::cx_mat> m_uSlices;              /**< Vector of Hamiltonian coupling matrices of slices*/
  std::vector<arma::cx_mat> m_suSlices;             /**< Vector of overlap coupling matrices of slices*/
  std::vector<std::array<double,2> > m_boundaries;  /**< Vector of arrays denoting boundaries*/
  double m_minSliceDist;                            /**< Minimum distance between slice boundaries*/
  arma::mat m_HMat;                                 /**< Device Hamiltonian */
  arma::mat m_SMat;                                 /**< Device overlap matrix */
  int m_llucAtoms;                                  /**< Number of atoms per unit cell in left lead */
  int m_rlucAtoms;                                  /**< Number of atoms per unit cell in right lead */
};
#endif
