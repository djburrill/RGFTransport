#include <iostream>
#include <string>
#include <readInput.h>
#include <deviceSlice.h>
#include <green.h>
#include <output.h>
#include <armadillo>
#include <cmath>
#include <vector>

/** Fermi-Dirac distribution.
*/
double fermi(double a_energy,double a_kT){
  return 1.0/(1.0+std::exp(a_energy/a_kT));
}

/** Main program function.
*/
int main(int argc, char* argv[]) {
  // Throw error if no file name specified
  if (argc == 1){
    std::cout << "Please specify an input file.";
    return 1;
  }

  // Variables
  std::string inputFileName = argv[1];
  std::vector<double> voltageVec;
  std::vector<double> currentVec;
  std::vector<double> DOSVec;
  std::vector<double> energyVec;
  std::vector<double> sgfVec_real;
  std::vector<double> sgfVec_imag;
  double charge = 1.602e-19;
  double planck = 6.626e-34;
  double coeff = 2*charge/planck;
  double kT = 0.1;
  output outputVals;

  //std::cout << "STARTED" << std::endl;

  // Parse input file
  readInput inputParams(inputFileName);

  //std::cout << "PARSED INPUT" << std::endl;

  // Slice device atoms
  deviceSlice dSlice(inputParams.getMolecule(),
                     inputParams.getDHMat(),
                     inputParams.getDSMat(),
                     inputParams.m_minSliceDist,
                     inputParams.m_llucAtoms,
                     inputParams.m_rlucAtoms);

  //std::cout << "SLICED DEVICE" << std::endl;

  // Set up Greens function class
  green GreensFuncs(inputParams.m_imag_eta,
                    inputParams.m_sc_conv,
                    inputParams.m_numBasisUCL,
                    inputParams.m_numBasisUCR,
                    inputParams.getLHMat(),
                    inputParams.getRHMat(),
                    inputParams.getDHMat(),
                    inputParams.getLSMat(),
                    inputParams.getRSMat(),
                    inputParams.getDSMat(),
                    dSlice.getHamSlices(),
                    dSlice.getOverSlices(),
                    dSlice.getUSlices(),
                    dSlice.getSUSlices());

  //std::cout << "SET UP GF CLASS" << std::endl;

  // Iterate over voltages
  //GreensFuncs.transmission(0.0);


  for (double volt=inputParams.m_vlow;volt<=inputParams.m_vhigh; volt+=inputParams.m_dv){
    // Variables
    arma::cx_double current = 0.0;
    std::vector<arma::cx_double> kernel;
    voltageVec.push_back(volt);

    // Iterate over energies
    for (double en=inputParams.m_elow; en<=inputParams.m_ehigh; en+=inputParams.m_de){
      // Compute transmission matrix
      GreensFuncs.transmission(en);

      // Take trace of transmission matrix
      arma::cx_double transmit = arma::trace(GreensFuncs.getTransMat());
      //std::cout << transmit << std::endl;

      // Compute Fermi functions
      double f1 = fermi(en-volt/2.0,kT);
      double f2 = fermi(en+volt/2.0,kT);
      arma::cx_double fDiff = f1-f2;
      //std::cout << "FDIFF: " << fDiff << std::endl;

      // Calculate kernel
      kernel.push_back(fDiff*transmit);

      // DEBUG
      //std::cout << sgf << std::endl;
      // Compute eigen-diagonal entries
      //arma::cx_mat sgL = GreensFuncs.getSGL();
      //std::cout << sgL << std::endl;
      //arma::cx_vec eigval = arma::eig_gen(sgL);
      //std::cout << eigval[0].real() << " " << eigval[0].imag() << " " << eigval[1].real() << " " << //eigval[1].imag() << " " << eigval[2].real() << " " << eigval[2].imag() << " " << //eigval[3].real() << " " << eigval[3].imag() << std::endl;
    }

    // Perform integration (trapezoidal)
    for (int i=0; i<kernel.size()-1; i++){
      current += 0.5*(kernel[i]+kernel[i+1])*inputParams.m_de;
    }

    // Save current
    currentVec.push_back(current.real());
    std::cout << volt << " " << current << std::endl;
  }


  /*
  // Calculate surface Green's function
  arma::cx_mat a_surfMat;
  arma::cx_mat a_uMat;
  arma::mat LHMat = inputParams.getLHMat();
  arma::mat LSMat = inputParams.getLSMat();
  for (double en=inputParams.m_elow; en<=inputParams.m_ehigh; en+=inputParams.m_de){
    GreensFuncs.surfGreen(a_surfMat,
                          a_uMat,
                          en,
                          inputParams.m_imag_eta,
                          inputParams.m_sc_conv,
                          inputParams.m_numBasisUCL,
                          LHMat,
                          LSMat);

    std::cout << "ENERGY: " << en << std::endl;
    std::cout << a_surfMat << std::endl;

    //arma::cx_vec eigval = arma::eig_gen(a_surfMat);

    // Compute DOS and save value
    //energyVec.push_back(en);
    //DOSVec.push_back(eigval[0].real());
  }
  */


  // Calculate DOS
  for (double en=inputParams.m_elow; en<=inputParams.m_ehigh; en+=inputParams.m_de){
    // Compute DOS and save value
    energyVec.push_back(en);
    DOSVec.push_back(GreensFuncs.calc_DOS(en));
  }

  // Calculate surface Green's function
  std::vector<arma::cx_mat> sgf = GreensFuncs.calc_SGF_EnRange(inputParams.m_elow,inputParams.m_ehigh,inputParams.m_de,0);

  // Calculate eigenvalues for surface Green's functions
  for (int index=0; index<sgf.size(); index++){
    arma::cx_vec eigval = arma::eig_gen(sgf[index]);
    sgfVec_real.push_back(eigval[0].real());
    sgfVec_imag.push_back(eigval[0].imag());
  }

  // Output IV curve
  outputVals.writeData("IV.dat",voltageVec,currentVec);

  // Output DOS
  outputVals.writeData("DOS.dat",energyVec,DOSVec);

  // Output sgf
  outputVals.writeData("sgf.dat",energyVec,sgfVec_real,sgfVec_imag);

}
