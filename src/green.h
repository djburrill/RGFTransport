/** Collection of methods for computing Green's functions.
*/
#ifndef _green_h_
#define _green_h_
#include <iostream>
#include <armadillo>
#include <complex>

class green{
public:
  /** Default constructor
  */
  green(){};

  /** Initialization constructor
  */
  green(double a_eta,
        double a_scConvergence,
        int a_basisPerUCL,
        int a_basisPerUCR,
        const arma::mat& a_HMatL,
        const arma::mat& a_HMatR,
        const arma::mat& a_HMatD,
        const arma::mat& a_SMatL,
        const arma::mat& a_SMatR,
        const arma::mat& a_SMatD,
        const std::vector<arma::cx_mat>& a_HamSlices,
        const std::vector<arma::cx_mat>& a_overSlices,
        const std::vector<arma::cx_mat>& a_uSlices,
        const std::vector<arma::cx_mat>& a_suSlices);

  /** Calculate left surface Green's function
  */
  void surfGreen_Left(double a_energy);

  /** Calculate right surface Green's function
  */
  void surfGreen_Right(double a_energy);

  /** Calculate surface Green's function
  */
  void surfGreen(arma::cx_mat& a_surfMat,
                 arma::cx_mat& a_uMat,
                 double a_energy,
                 double a_eta,
                 double a_scConvergence,
                 int a_basisPerUC,
                 arma::mat& a_HMat,
                 arma::mat& a_SMat);

  /** Calculate left self energy matrix
  */
  void selfEnergy_Left();

  /** Calculate right self energy matrix
  */
  void selfEnergy_Right();

  /** Calculate self energy matrix
  */
  void selfEnergy(arma::cx_mat& a_selfEnergyMat,
                  arma::cx_mat& a_surfMat,
                  arma::cx_mat& a_uMat);

  /** Calculate left broadening matrix
  */
  void broadening_Left();

  /** Calculate right broadening matrix
  */
  void broadening_Right();

  /** Calculate broadening matrix
  */
  void broadening(arma::cx_mat& a_broadeningMat,
                  arma::cx_mat& a_selfEnergyMat);

  /** Get left surface Green's function
  */
  arma::cx_mat getSGL(){return m_sgL;};

  /** Get right surface Green's function
  */
  arma::cx_mat getSGR(){return m_sgR;};

  /** Calculate device Green's function matrices
  */
  void deviceGF(arma::cx_mat& a_DeviceGF_Left,
                arma::cx_mat& a_DeviceGF_Right,
                double a_energy,
                std::vector<arma::cx_mat>& a_HamSlices,
                std::vector<arma::cx_mat>& a_overSlices,
                std::vector<arma::cx_mat>& a_uSlices,
                std::vector<arma::cx_mat>& a_suSlices,
                arma::cx_mat& a_uL,
                arma::cx_mat& a_uR,
                arma::cx_mat& a_sgL,
                arma::cx_mat& a_sgR);

  /** Calculate device Green's function matrices - INTERNAL
  */
  void deviceGF(double a_energy);

  /** Calculate transmission matrix - INTERNAL
  */
  void transmission(double a_energy);

  /** Get transmission matrix
  */
  arma::cx_mat getTransMat(){return m_TransMat;};

  /** Calculate density of states - INTERNAL
  */
  double calc_DOS(double a_energy);

  /** Calculate surface Green's function over an energy range
  */
  std::vector<arma::cx_mat> calc_SGF_EnRange(double a_enStart, double a_enStop, double a_enStep, int a_side);

private:
  // Private variables
  double m_eta;             /**< Imaginary energy factor */
  double m_scConvergence;   /**< Self-consisten convergence paramter */
  int m_basisPerUCL;        /**< Number of basis function per unit cell left lead */
  int m_basisPerUCR;        /**< Number of basis function per unit cell right lead */
  arma::mat m_HMatL;        /**< Left lead Hamiltonian matrix */
  arma::mat m_HMatR;        /**< Right lead Hamiltonian matrix */
  arma::mat m_HMatD;        /**< Device Hamiltonian matrix */
  arma::mat m_SMatL;        /**< Left lead overlap matrix */
  arma::mat m_SMatR;        /**< Right lead overlap matrix */
  arma::mat m_SMatD;        /**< Device overlap matrix */
  arma::cx_mat m_sgL;       /**< Left surface Greens function */
  arma::cx_mat m_sgR;       /**< Right surface Greens function */
  arma::cx_mat m_uL;        /**< Left coupling matrix */
  arma::cx_mat m_uR;        /**< Right coupling matrix */
  arma::cx_mat m_seL;       /**< Left self energy matrix */
  arma::cx_mat m_seR;       /**< Right self energy matrix */
  arma::cx_mat m_bL;        /**< Left broadening matrix */
  arma::cx_mat m_bR;        /**< Right broadening matrix */
  std::vector<arma::cx_mat> m_HamSlices;     /**< Vector of Hamiltonians of slices */
  std::vector<arma::cx_mat> m_overSlices;    /**< Vector of overlap matrices of slices */
  std::vector<arma::cx_mat> m_uSlices;       /**< Vector of Hamiltonian coupling matrices of slices */
  std::vector<arma::cx_mat> m_suSlices;      /**< Vector of overlap coupling matrices of slices */
  arma::cx_mat m_DeviceGF_Left;              /**< Device Green's function - Left */
  arma::cx_mat m_DeviceGF_Right;             /**< Device Green's function - Right */
  arma::cx_mat m_TransMat;                   /**< Transmission matrix */
};
#endif
