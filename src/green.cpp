#include <green.h>

green::green(double a_eta,
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
             const std::vector<arma::cx_mat>& a_suSlices){
  // Set variables
  m_eta = a_eta;
  m_scConvergence = a_scConvergence;
  m_basisPerUCL = a_basisPerUCL;
  m_basisPerUCR = a_basisPerUCR;
  m_HMatL = a_HMatL;
  m_HMatR = a_HMatR;
  m_HMatD = a_HMatD;
  m_SMatL = a_SMatL;
  m_SMatR = a_SMatR;
  m_SMatD = a_SMatD;
  m_HamSlices = a_HamSlices;
  m_overSlices = a_overSlices;
  m_uSlices = a_uSlices;
  m_suSlices = a_suSlices;
}

void green::surfGreen_Left(double a_energy){
  surfGreen(m_sgL,m_uL,a_energy,m_eta,m_scConvergence,m_basisPerUCL,m_HMatL,m_SMatL);
}

void green::surfGreen_Right(double a_energy){
  //surfGreen(m_sgR,m_uR,a_energy,m_eta,m_scConvergence,m_basisPerUCR,m_HMatR,m_SMatR);
  //m_uR = arma::trans(m_uR);
  m_sgR = m_sgL;
  m_uR = arma::trans(m_uL);
}

void green::surfGreen(arma::cx_mat& a_surfMat,
                      arma::cx_mat& a_uMat,
                      double a_energy,
                      double a_eta,
                      double a_scConvergence,
                      int a_basisPerUC,
                      arma::mat& a_HMat,
                      arma::mat& a_SMat){
  // Variables
  int numDivs = arma::size(a_HMat)[0]/a_basisPerUC;
  int hCell = numDivs/2.0;
  int colStart = hCell*(a_basisPerUC);
  int colEnd = ((hCell+1)*a_basisPerUC)-1;
  int ucolStart = (hCell+1)*(a_basisPerUC);
  int ucolEnd = ((hCell+2)*a_basisPerUC)-1;
  int maxCycles = 100;
  int numCycles = 0;
  arma::cx_double cmplxEnergy(a_energy,a_eta);
  arma::cx_mat HMat = arma::conv_to<arma::cx_mat>::from(a_HMat);
  arma::cx_mat SMat = arma::conv_to<arma::cx_mat>::from(a_SMat);

  // Identify h and u matrices
  arma::cx_mat hMat = HMat(arma::span(colStart,colEnd),arma::span(colStart,colEnd));
  arma::cx_mat sMat = SMat(arma::span(colStart,colEnd),arma::span(colStart,colEnd));
  arma::cx_mat uMat = HMat(arma::span(colStart,colEnd),arma::span(ucolStart,ucolEnd));
  arma::cx_mat suMat = SMat(arma::span(colStart,colEnd),arma::span(ucolStart,ucolEnd));
  a_uMat = uMat;

  //std::cout << "HAMILTONIAN: " << hMat << std::endl;
  //std::cout << "COUPLING: " << uMat << std::endl;

  // Calculate initial components
  //arma::cx_mat invEn = arma::inv(cmplxEnergy*sMat - hMat);
  arma::cx_mat invEn = arma::inv(cmplxEnergy*arma::eye(arma::size(hMat)) - hMat);
  //arma::cx_mat alphaOld = (cmplxEnergy*suMat - uMat)*invEn*(cmplxEnergy*suMat - uMat);
  //arma::cx_mat betaOld = arma::trans(cmplxEnergy*suMat - uMat)*invEn*arma::trans(cmplxEnergy*suMat - uMat);
  //arma::cx_mat epSOld = hMat + (cmplxEnergy*suMat - uMat)*invEn*arma::trans(cmplxEnergy*suMat - uMat);
  //arma::cx_mat epOld = hMat + (cmplxEnergy*suMat - uMat)*invEn*arma::trans(cmplxEnergy*suMat - uMat) + arma::trans(cmplxEnergy*suMat - uMat)*invEn*(cmplxEnergy*suMat - uMat);
  //arma::cx_mat alphaOld = uMat;
  //arma::cx_mat betaOld = arma::trans(uMat);
  //arma::cx_mat epSOld = hMat;
  //arma::cx_mat epOld = epSOld;
  //arma::cx_mat alphaOld = -(cmplxEnergy*suMat - uMat);
  //arma::cx_mat betaOld = -arma::trans(cmplxEnergy*suMat - uMat);
  //arma::cx_mat epSOld = -(cmplxEnergy*sMat - hMat);
  //arma::cx_mat epOld = epSOld;
  arma::cx_mat alphaOld = uMat*invEn*uMat;
  arma::cx_mat betaOld = uMat*invEn*arma::trans(uMat);
  arma::cx_mat epSOld = hMat + uMat*invEn*arma::trans(uMat);
  arma::cx_mat epOld = epSOld + arma::trans(uMat)*invEn*uMat;

  // Dummy updates
  arma::cx_mat alphaNew = 3*a_scConvergence*alphaOld;
  arma::cx_mat betaNew = 3*a_scConvergence*betaOld;
  arma::cx_mat epSNew = epOld + alphaOld*invEn*betaOld;
  arma::cx_mat epNew = epSNew + betaOld*invEn*alphaOld;
  double alphaDiff = arma::norm(alphaNew-alphaOld);
  double betaDiff = arma::norm(betaNew-betaOld);
  double epSDiff = arma::norm(epSNew-epSOld);
  double epDiff = arma::norm(epNew-epOld);

  // Self-consistent cycle
  while (((alphaDiff>a_scConvergence) || (betaDiff>a_scConvergence) || (epSDiff>a_scConvergence) || (epDiff>a_scConvergence)) && (numCycles<maxCycles)){
    //std::cout << "ALPHA OLD: " << alphaOld << std::endl;
    //std::cout << "BETA OLD: " << betaOld << std::endl;
    //std::cout << "EPS OLD: " << epSOld << std::endl;
    //std::cout << "EP OLD: " << epOld << std::endl;

    // Calculate updates
    //invEn = arma::inv(cmplxEnergy*sMat - epOld);
    invEn = arma::inv(cmplxEnergy*arma::eye(arma::size(epOld)) - epOld);
    alphaNew = alphaOld*invEn*alphaOld;
    betaNew = betaOld*invEn*betaOld;
    epNew = epOld + betaOld*invEn*alphaOld + alphaOld*invEn*betaOld;
    epSNew = epSOld + alphaOld*invEn*betaOld;
    //epSNew = epOld + alphaOld*invEn*betaOld;
    //epNew = epSNew + betaOld*invEn*alphaOld;

    // Calculate differences
    alphaDiff = arma::norm(alphaNew-alphaOld);
    betaDiff = arma::norm(betaNew-betaOld);
    epSDiff = arma::norm(epSNew-epSOld);
    epDiff = arma::norm(epNew-epOld);

    //std::cout << "ALPHA NEW: " << alphaNew << std::endl;
    //std::cout << "BETA NEW: " << betaNew << std::endl;
    //std::cout << "EPS NEW: " << epSNew << std::endl;
    //std::cout << "EP NEW: " << epNew << std::endl;
//
    //std::cout << "ONSITE g\n" << invEn << std::endl;
    //std::cout << "ALPHA DIFF: " << alphaDiff << std::endl;
    //std::cout << "BETA DIFF: " << betaDiff << std::endl;
    //std::cout << "EPS DIFF: " << epSDiff << std::endl;
    //std::cout << "EP DIFF: " << epDiff << std::endl;

    // Update old terms
    alphaOld = alphaNew;
    betaOld = betaNew;
    epSOld = epSNew;
    epOld = epNew;

    // Update cycle number
    numCycles += 1;
  }

  // Catch max cycles
  if (numCycles >= maxCycles){
    std::cout << "MAX SELF CONSISTENT CYCLES REACHED" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Calculate surface Green's function
  //a_surfMat = epSNew;
  //a_surfMat = arma::inv(cmplxEnergy*sMat - epSNew);
  //std::cout << "SURFACE: " << epSNew << std::endl;
  //std::cout << "BULK: " << epNew << std::endl;
  a_surfMat = arma::inv(cmplxEnergy*arma::eye(arma::size(epSNew)) - epSNew);
}

void green::selfEnergy_Left(){
  // Calculate left self energy
  selfEnergy(m_seL,m_sgL,m_uL);
}

void green::selfEnergy_Right(){
  // Calculate left self energy
  arma::cx_mat tmpUR = arma::trans(m_uR);
  selfEnergy(m_seR,m_sgR,tmpUR);
}

void green::selfEnergy(arma::cx_mat& a_selfEnergyMat,
                       arma::cx_mat& a_surfMat,
                       arma::cx_mat& a_uMat){
  // Calcualte self energy
  //std::cout << a_uMat << std::endl;
  //std::cout << a_surfMat << std::endl;

  a_selfEnergyMat = a_uMat*a_surfMat*arma::trans(a_uMat);
}

void green::broadening_Left(){
  broadening(m_bL,m_seL);
}

void green::broadening_Right(){
  broadening(m_bR,m_seR);
}

void green::broadening(arma::cx_mat& a_broadeningMat,
                       arma::cx_mat& a_selfEnergyMat){
  // Calculate broadening matrix
  arma::cx_double unitImag(0.0,1.0);

  //std::cout << "SELF ENERGY\n" << a_selfEnergyMat << std::endl;
  //std::cout << "SELF ENERGY TRANSPOSE\n" << arma::trans(a_selfEnergyMat) << std::endl;
  a_broadeningMat = unitImag*(a_selfEnergyMat-arma::trans(a_selfEnergyMat));
}

void green::deviceGF(arma::cx_mat& a_DeviceGF_Left,
                     arma::cx_mat& a_DeviceGF_Right,
                     double a_energy,
                     std::vector<arma::cx_mat>& a_HamSlices,
                     std::vector<arma::cx_mat>& a_overSlices,
                     std::vector<arma::cx_mat>& a_uSlices,
                     std::vector<arma::cx_mat>& a_suSlices,
                     arma::cx_mat& a_uL,
                     arma::cx_mat& a_uR,
                     arma::cx_mat& a_sgL,
                     arma::cx_mat& a_sgR){
  // Variables
  int numIters = a_HamSlices.size() + 1;    /**< Number of iterations for calculating Green's functions */
  arma::cx_mat tmpGreen_Left;
  arma::cx_mat tmpGreen_Right;
  arma::cx_mat tmpGreen_OnSite;

  /*
  //DEBUG
  std::cout << "SLICES" << std::endl;
  for (int index=0; index<a_uSlices.size(); index++){
    std::cout << a_uSlices[index] << std::endl;
  }
  */

  // Left Green's function
  for (int iteration=0; iteration<numIters; iteration++){
    //std::cout << iteration << std::endl;
    // First iteration
    if (iteration == 0){

      /*
      // DEBUG
      std::cout << "LEFT OVERLAP: " << iteration << " " << a_overSlices[iteration] << std::endl;
      std::cout << "LEFT HAMILTONIAN: " << iteration << " " << a_HamSlices[iteration] << std::endl;
      std::cout << "LEFT COUPLING: " << iteration << " " << a_uL << std::endl;
      std::cout << "LEFT CGF: " << iteration << " " << arma::trans(a_uL)*a_sgL*a_uL << std::endl;

      std::cout << "LEFT TRANS COUPLING: " << iteration << " " <<  arma::size(arma::trans(a_uL)) << std::endl;
      std::cout << "LEFT SGL: " << iteration << " " <<  arma::size(a_sgL) << std::endl;
      std::cout << "LEFT COUPLING: " << iteration << " " <<  arma::size(a_uL) << std::endl;
      */

      // Compute on-site Green's function matrix
      tmpGreen_OnSite = a_overSlices[iteration]*a_energy - a_HamSlices[iteration] - arma::trans(a_uL)*a_sgL*a_uL;
      tmpGreen_OnSite = arma::inv(tmpGreen_OnSite);

      // Compute iterative Green's function matrix
      tmpGreen_Left = a_sgL*a_uL*tmpGreen_OnSite;
      //std::cout << "LEFT GF: " << iteration << " " <<  tmpGreen_Left << std::endl;

      // Continue to next iteration
      //continue;
    }
    // Last iteration
    else if (iteration == numIters-1){
      // Compute iterative Green's function matrix
      tmpGreen_Left = tmpGreen_Left*a_uR*a_sgR;

      // DEBUG
      //std::cout << "LEFT GF: " << iteration << " " << tmpGreen_Left << std::endl;
      //std::cout << "LEFT SGF: " << iteration << " " << a_sgL << std::endl;

      // Continue to next iteration
      //continue;
    }
    // All other iterations
    else{

      /*
      // DEBUG
      std::cout << "LEFT OVERLAP: " << iteration << " " << a_overSlices[iteration] << std::endl;
      std::cout << "LEFT HAMILTONIAN: " << iteration << " " << a_HamSlices[iteration] << std::endl;
      std::cout << "LEFT COUPLING: " << iteration << " " << a_uSlices[iteration-1]<< std::endl;
      std::cout << "LEFT CGF: " << iteration << " " << arma::trans(a_uSlices[iteration-1])*tmpGreen_OnSite*a_uSlices[iteration-1] << std::endl;

      std::cout << "LEFT TRANS COUPLING: " << iteration << " " <<  arma::size(arma::trans(a_uSlices[iteration-1])) << std::endl;
      std::cout << "LEFT ONSITE: " << iteration << " " <<  arma::size(tmpGreen_OnSite) << std::endl;
      std::cout << "LEFT COUPLING: " << iteration << " " <<  arma::size(a_uSlices[iteration-1]) << std::endl;
      */

      // Compute on-site Green's function matrix
      tmpGreen_OnSite = a_overSlices[iteration]*a_energy - a_HamSlices[iteration] - arma::trans(a_uSlices[iteration-1])*tmpGreen_OnSite*a_uSlices[iteration-1];
      tmpGreen_OnSite = arma::inv(tmpGreen_OnSite);

      // Compute iterative Green's function matrix
      tmpGreen_Left = tmpGreen_Left*a_uSlices[iteration-1]*tmpGreen_OnSite;
      //std::cout << "LEFT GF: " << iteration << " " <<  tmpGreen_Left << std::endl;

      // Continue to next iteration
      //continue;
    }

    /*
    // DEBUG
    std::cout << "LEFT ONSITE: " << iteration << " " << tmpGreen_OnSite << std::endl;
    std::cout << "LEFT GF ITERATION: " << iteration << " " << tmpGreen_Left << std::endl;
    */
  }

  // Right Green's function
  for (int iteration=0; iteration<numIters; iteration++){
    // Easy tracking index
    int sliceIndex = numIters-2-iteration;

    // First iteration
    if (iteration == 0){

      /*
      // DEBUG
      std::cout << "RIGHT OVERLAP: " << iteration << " " << a_overSlices[sliceIndex] << std::endl;
      std::cout << "RIGHT HAMILTONIAN: " << iteration << " " << a_HamSlices[sliceIndex] << std::endl;
      std::cout << "RIGHT COUPLING: " << iteration << " " << a_uR << std::endl;
      std::cout << "RIGHT CGF: " << iteration << " " << arma::trans(a_uR)*a_sgR*a_uR << std::endl;

      std::cout << "RIGHT TRANS COUPLING: " << iteration << " " <<  arma::size(arma::trans(a_uR)) << std::endl;
      std::cout << "RIGHT SGL: " << iteration << " " <<  arma::size(a_sgR) << std::endl;
      std::cout << "RIGHT COUPLING: " << iteration << " " <<  arma::size(a_uR) << std::endl;
      */

      // Compute on-site Green's function matrix
      tmpGreen_OnSite = a_overSlices[sliceIndex]*a_energy - a_HamSlices[sliceIndex] - arma::trans(a_uR)*a_sgR*a_uR;
      //tmpGreen_OnSite = a_overSlices[sliceIndex]*a_energy - a_HamSlices[sliceIndex] - a_uR*a_sgR*arma::trans(a_uR);
      tmpGreen_OnSite = arma::inv(tmpGreen_OnSite);

      // Compute iterative Green's function matrix
      tmpGreen_Right = a_sgR*a_uR*tmpGreen_OnSite;
      //std::cout << "RIGHT GF: " << iteration << " " <<  tmpGreen_Right << std::endl;
      //tmpGreen_Right = a_sgR*arma::trans(a_uR)*tmpGreen_OnSite;

      // Continue to next iteration
      //continue;
    }
    // Last iteration
    else if (iteration == numIters-1){
      // Compute iterative Green's function matrix
      tmpGreen_Right = tmpGreen_Right*a_uR*a_sgL;
      //tmpGreen_Right = tmpGreen_Right*arma::trans(a_uR)*a_sgL;

      // DEBUG
      //std::cout << "RIGHT GF: " << iteration << " " << tmpGreen_Right << std::endl;
      //std::cout << "RIGHT SGF: " << iteration << " " << a_sgR << std::endl;

      // Continue to next iteration
      //continue;
    }
    // All other iterations
    else{

      /*
      // DEBUG
      std::cout << "RIGHT OVERLAP: " << iteration << " " << a_overSlices[sliceIndex] << std::endl;
      std::cout << "RIGHT HAMILTONIAN: " << iteration << " " << a_HamSlices[sliceIndex] << std::endl;
      std::cout << "SLICE INDEX: " << sliceIndex << std::endl;
      std::cout << "RIGHT COUPLING: " << iteration << " " << a_uSlices[sliceIndex] << std::endl;
      std::cout << "RIGHT CGF: " << iteration << " " << arma::trans(a_uSlices[sliceIndex])*tmpGreen_OnSite*a_uSlices[sliceIndex] << std::endl;

      std::cout << "RIGHT TRANS COUPLING: " << iteration << " " <<  arma::size(a_uSlices[sliceIndex]) << std::endl;
      std::cout << "RIGHT ONSITE: " << iteration << " " <<  arma::size(tmpGreen_OnSite) << std::endl;
      std::cout << "RIGHT COUPLING: " << iteration << " " <<  arma::size(arma::trans(a_uSlices[sliceIndex])) << std::endl;
      */

      //std::cout << "SLICE INDEX: " << sliceIndex << std::endl;


      // Compute on-site Green's function matrix
      //tmpGreen_OnSite = a_overSlices[sliceIndex]*a_energy - a_HamSlices[sliceIndex] - arma::trans(a_uSlices[sliceIndex])*tmpGreen_OnSite*a_uSlices[sliceIndex];
      //tmpGreen_OnSite = a_overSlices[sliceIndex]*a_energy - a_HamSlices[sliceIndex] - arma::trans(a_uSlices[sliceIndex])*tmpGreen_OnSite*a_uSlices[sliceIndex];
      tmpGreen_OnSite = a_overSlices[sliceIndex]*a_energy - a_HamSlices[sliceIndex] - a_uSlices[sliceIndex]*tmpGreen_OnSite*arma::trans(a_uSlices[sliceIndex]);
      tmpGreen_OnSite = arma::inv(tmpGreen_OnSite);

      // Compute iterative Green's function matrix
      //tmpGreen_Right = tmpGreen_Right*a_uSlices[sliceIndex-1]*tmpGreen_OnSite;
      //tmpGreen_Right = tmpGreen_Right*a_uSlices[sliceIndex]*tmpGreen_OnSite;
      tmpGreen_Right = tmpGreen_Right*arma::trans(a_uSlices[sliceIndex])*tmpGreen_OnSite;
      //std::cout << "RIGHT GF: " << iteration << " " <<  tmpGreen_Right << std::endl;

      // Continue to next iteration
      //continue;
    }

    /*
    // DEBUG
    std::cout << "RIGHT ONSITE: " << iteration << " " << tmpGreen_OnSite << std::endl;
    std::cout << "RIGHT GF ITERATION: " << iteration << " " << tmpGreen_Right << std::endl;
    */
  }

  // Write Green's function matrices
  a_DeviceGF_Left = tmpGreen_Left;
  a_DeviceGF_Right = tmpGreen_Right;
}

void green::deviceGF(double a_energy){
  deviceGF(m_DeviceGF_Left,m_DeviceGF_Right,a_energy,m_HamSlices,m_overSlices,m_uSlices,m_suSlices,m_uL,m_uR,m_sgL,m_sgR);

  /*
  std::cout << "LEFT COUPLING" << std::endl;
  std::cout << m_uL << std::endl;
  std::cout << "LEFT SGF" << std::endl;
  std::cout << m_sgL << std::endl;
  std::cout << "LEFT DEVICE" << std::endl;
  std::cout << m_DeviceGF_Left << std::endl;

  std::cout << "RIGHT COUPLING" << std::endl;
  std::cout << m_uR << std::endl;
  std::cout << "RIGHT SGF" << std::endl;
  std::cout << m_sgR << std::endl;
  std::cout << "RIGHT DEVICE" << std::endl;
  std::cout << m_DeviceGF_Right << std::endl;
  */
}

void green::transmission(double a_energy){
  // Variables

  // Calculate surface Green's functions
  surfGreen_Left(a_energy);
  surfGreen_Right(a_energy);
  //std::cout << "SURFACE GF CALCULATED" << std::endl;

  // Calculate self energy matrices
  selfEnergy_Left();
  selfEnergy_Right();
  //std::cout << "SE CALCULATED" << std::endl;

  //std::cout << m_seL << std::endl;
  //std::cout << m_seR << std::endl;

  // Calculate broadening matrices
  broadening_Left();
  broadening_Right();
  //std::cout << "BROADENING MATRICES CALCULATED" << std::endl;

  // Calculate device Green's function matrices
  deviceGF(a_energy);
  //std::cout << "DEVICE GF CALCULATED" << std::endl;


  //DEBUG
  //std::cout << "GAMMA LEFT\n" << m_bL << std::endl;
  //std::cout << "DEVICE LEFT\n" << m_DeviceGF_Left << std::endl;
  //std::cout << "GAMMA RIGHT\n " << m_bR << std::endl;
  //std::cout << "DEVICE RIGHT\n" << m_DeviceGF_Right << std::endl;


  // Calculate transmission matrix
  m_TransMat = m_bL*m_DeviceGF_Left*m_bR*arma::trans(m_DeviceGF_Right);
  //std::cout << "TRANSMISSION CALCULATED" << std::endl;
}

double green::calc_DOS(double a_energy){
  // Calculate surface Green's functions
  surfGreen_Left(a_energy);
  surfGreen_Right(a_energy);

  // Calculate device Green's function matrices
  deviceGF(a_energy);

  // Return DOS
  return (-1.0/M_PI)*arma::trace(m_DeviceGF_Left).imag();
}

std::vector<arma::cx_mat> green::calc_SGF_EnRange(double a_enStart, double a_enStop, double a_enStep, int a_side){
  // Variables
  std::vector<arma::cx_mat> outGF;

  // Loop over energies
  for (double en=a_enStart; en<=a_enStop; en+=a_enStep){
    // Initialize surface Greens function
    arma::cx_mat sgf;

    // Left side
    if (a_side == 0){
      surfGreen(sgf,m_uL,en,m_eta,m_scConvergence,m_basisPerUCL,m_HMatL,m_SMatL);
    }
    // Right side
    else if (a_side == 1){
      surfGreen(sgf,m_uR,en,m_eta,m_scConvergence,m_basisPerUCR,m_HMatR,m_SMatR);
    }

    // Append surface Green's function
    outGF.push_back(sgf);
  }

  return outGF;
}
