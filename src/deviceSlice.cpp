#include <deviceSlice.h>

deviceSlice::deviceSlice(const std::vector<atom>& a_molecule,
                         const arma::mat& a_HMat,
                         const arma::mat& a_SMat,
                         double a_minSliceDist,
                         int a_llucAtoms,
                         int a_rlucAtoms){
  // Set member variables
  m_molecule = a_molecule;
  m_minSliceDist = a_minSliceDist;
  m_HMat = a_HMat;
  m_SMat = a_SMat;
  m_llucAtoms = a_llucAtoms;
  m_rlucAtoms = a_rlucAtoms;

  // Slice molecule
  slice();

  // Slice device matrices
  //std::cout << "SLICING DEVICE MATRICES" << std::endl;
  sliceHam(m_HamSlices,m_overSlices,m_uSlices,m_suSlices,m_slices,m_HMat,m_SMat);

  /*
  std::cout << "SLICE HAMILTONIANS" << std::endl;
  for (int i=0; i<m_HamSlices.size(); i++){
    std::cout << "SLICE #" << i << std::endl;
    std::cout << m_HamSlices[i] << std::endl;
  }

  std::cout << "SLICE OVERLAPS" << std::endl;
  for (int i=0; i<m_overSlices.size(); i++){
    std::cout << "SLICE #" << i << std::endl;
    std::cout << m_overSlices[i] << std::endl;
  }

  std::cout << "SLICE HAMILTONIAN COUPLING" << std::endl;
  for (int i=0; i<m_uSlices.size(); i++){
    std::cout << "SLICE #" << i << "-" << i+1 << std::endl;
    std::cout << m_uSlices[i] << std::endl;
  }

  std::cout << "SLICE OVERLAP COUPLING" << std::endl;
  for (int i=0; i<m_suSlices.size(); i++){
    std::cout << "SLICE #" << i << "-" << i+1 << std::endl;
    std::cout << m_suSlices[i] << std::endl;
  }
  */

  /*
  for (int i=0; i<m_slices.size(); i++){
    for (int j=0; j<m_slices[i].size(); j++){
      //std::cout << m_slices[i][j].getAtomNum() << std::endl;
      //std::cout << m_slices[i][j].getPos()[0] << std::endl;
    }
  }
  */
}

void deviceSlice::slice(){
  // Variables
  std::array<double,2> boundEnds;
  int numDivs;
  double eta = 0.01;                /**< Slightly increases the size of slices to avoid problems with checking 'double' equalities*/
  std::vector<atom> leftLead;
  std::vector<atom> rightLead;

  // Remove lead unit cell atoms from the ends of the device
  // Left lead
  for (int i=0; i<m_llucAtoms; i++){
    // Add lead atom
    leftLead.push_back(m_molecule[0]);

    // Remove lead from molecule
    m_molecule.erase(m_molecule.begin());
  }

  // Right lead
  for (int i=0; i<m_llucAtoms; i++){
    // Add lead atom
    int molSize = m_molecule.size();
    rightLead.push_back(m_molecule[molSize-1]);

    // Remove lead from molecule
    m_molecule.erase(m_molecule.end()-1);
  }

  // Determine global device boundaries
  for (int i=0; i<m_molecule.size(); i++){
    double xPos = m_molecule[i].getPos()[0];

    // First iteration set boundaries
    if (i == 0){
      boundEnds[0] = xPos;
      boundEnds[1] = xPos;
      continue;
    }

    // Check low end
    if (xPos < boundEnds[0]){
      boundEnds[0] = xPos;
      continue;
    }

    // Check high end
    if (xPos > boundEnds[1]){
      boundEnds[1] = xPos;
      continue;
    }
  }

  // Determine number of divisions
  numDivs = (int) std::floor(((boundEnds[1]-boundEnds[0])/m_minSliceDist));

  if (numDivs < 1){
    numDivs = 0;
  }

  // Set up slice boundaries
  for (int i=0; i<=numDivs; i++){
    // Lower bound
    double lowerBound = (boundEnds[0]-eta)+i*m_minSliceDist;
    double upperBound;

    // Upper bound
    if (i+1 > numDivs){
      upperBound = boundEnds[1]+eta;
    }
    else{
      upperBound = lowerBound + m_minSliceDist;
    }

    // Add boundary set
    std::array<double,2> tmpArray{lowerBound,upperBound};
    m_boundaries.push_back(tmpArray);
  }

  // Set size of m_slices
  if (m_molecule.size() == 0){
    // Adds leads to slices
    m_slices.push_back(rightLead);
    m_slices.insert(m_slices.begin(),leftLead);
  }
  else{
    m_slices.resize(m_boundaries.size());

    // Add atoms to slices
    for (int i=0; i<m_molecule.size(); i++){
      int slicePos = (int) std::floor((m_molecule[i].getPos()[0]-boundEnds[0])/m_minSliceDist);

      m_slices[slicePos].push_back(m_molecule[i]);
    }

    // Remove empty slices
    for (int i=0; i<m_slices.size()-1; i++){
      int sizeNextSlice = m_slices[i+1].size();   /**< Number of atoms in next slice*/

      // If next slice has no atoms
      if (sizeNextSlice == 0){
        // Set correct boundaries
        m_boundaries[i][1] = m_boundaries[i+1][1];

        // Remove empty elements
        m_slices.erase(m_slices.begin()+i+1);
        m_boundaries.erase(m_boundaries.begin()+i+1);

        // Set i back one step
        i--;
      }
    }

    // Adds leads to slices
    m_slices.push_back(rightLead);
    m_slices.insert(m_slices.begin(),leftLead);
  }
}

void deviceSlice::sliceHam(std::vector<arma::cx_mat>& a_HamSlices,
                           std::vector<arma::cx_mat>& a_overSlices,
                           std::vector<arma::cx_mat>& a_uSlices,
                           std::vector<arma::cx_mat>& a_suSlices,
                           std::vector<std::vector<atom> >& a_slices,
                           arma::mat& a_HMat,
                           arma::mat& a_SMat){
  // Variables
  int numSlices = a_slices.size();

  //std::cout << "HAMILTONIAN\n" << a_HMat << std::endl;
  //std::cout << "OVERLAP\n" << a_SMat << std::endl;
  //std::cout << "HAMILTONIAN SIZE: " << arma::size(a_HMat) << std::endl;
  //std::cout << "OVERLAP SIZE: " << arma::size(a_SMat) << std::endl;

  // Iterate through atom slices
  for (int sliceNum=0; sliceNum<numSlices; sliceNum++){
    // Determine size of temp matrix
    int tmpSize = 0;
    int numAtoms = a_slices[sliceNum].size();
    for (int i=0; i<numAtoms; i++){
      tmpSize += a_slices[sliceNum][i].getNumBasis();
    }

    // Set up temp slice matrices
    arma::cx_mat tmpHMat(tmpSize,tmpSize);
    arma::cx_mat tmpSMat(tmpSize,tmpSize);

    //std::cout << "TMPSIZE: " << tmpSize << std::endl;

    // Counter variables
    int offsetRow = 0;
    int offsetColumn = 0;

    // Iterate through atoms in slice
    for (int atom1=0; atom1<numAtoms; atom1++){
      // Get atom 1 index
      int atom1Idx = a_slices[sliceNum][atom1].getMatIdx();
      int atom1NumBasis = a_slices[sliceNum][atom1].getNumBasis();

      //std::cout << "ATOM1IDX: " << atom1Idx << std::endl;
      //std::cout << "ATOM1NUMBASIS: " << atom1NumBasis << std::endl;

      // Reset column offset
      offsetColumn = 0;

      // Iterate through secondary atoms
      for (int atom2=atom1; atom2<numAtoms; atom2++){
        // Get atom 2 index
        int atom2Idx = a_slices[sliceNum][atom2].getMatIdx();
        int atom2NumBasis = a_slices[sliceNum][atom2].getNumBasis();

        //std::cout << "ATOM2IDX: " << atom2Idx << std::endl;
        //std::cout << "ATOM2NUMBASIS: " << atom2NumBasis << std::endl;

        // Loop over atom 1 basis functions
        for (int a1Basis=0; a1Basis<atom1NumBasis; a1Basis++){
          // Loop over atom 2 basis functions
          for (int a2Basis=0; a2Basis<atom2NumBasis; a2Basis++){
            //std::cout << "TMPROW: " << offsetRow+a1Basis << std::endl;
            //std::cout << "TMPCOL: " << offsetColumn+a2Basis << std::endl;
            //std::cout << "ROW: " << atom1Idx+a1Basis << std::endl;
            //std::cout << "COL: " << atom2Idx+a2Basis << std::endl;

            // Set slice matrix elements
            tmpHMat(offsetRow+a1Basis,offsetColumn+a2Basis) = a_HMat(atom1Idx+a1Basis,atom2Idx+a2Basis);
            tmpSMat(offsetRow+a1Basis,offsetColumn+a2Basis) = a_SMat(atom1Idx+a1Basis,atom2Idx+a2Basis);

            //std::cout << "SET MATRIX ELEMENTS" << std::endl;
          }
        }

        // Increment column
        offsetColumn += atom2NumBasis;
      }

      // Increment row
      offsetRow += atom1NumBasis;
    }

    // Add matrix slices
    a_HamSlices.push_back(tmpHMat);
    a_overSlices.push_back(tmpSMat);
    //std::cout << "SETTING UP COUPLING MATRICES" << std::endl;

    // Set up coupling matrices
    // Only run if not the last slice
    if (sliceNum < numSlices-1){
      // Determine size of temp matrices
      int numAtoms1 = a_slices[sliceNum].size();
      int numAtoms2 = a_slices[sliceNum+1].size();
      int tmpSizeRow = 0;
      int tmpSizeColumn = 0;

      for (int i=0; i<numAtoms1; i++){
        tmpSizeRow += a_slices[sliceNum][i].getNumBasis();
      }

      for (int i=0; i<numAtoms2; i++){
        tmpSizeColumn += a_slices[sliceNum+1][i].getNumBasis();
      }

      // Set up matrices
      arma::cx_mat tmpUMat(tmpSizeRow,tmpSizeColumn);
      arma::cx_mat tmpSUMat(tmpSizeRow,tmpSizeColumn);

      //std::cout << "TMPSIZE ROW: " << tmpSizeRow << std::endl;
      //std::cout << "TMPSIZE COLUMN: " << tmpSizeColumn << std::endl;

      // Counter variables
      int offsetRow = 0;
      int offsetColumn = 0;

      // Iterate through atoms in slice
      for (int atom1=0; atom1<numAtoms1; atom1++){
        // Get atom 1 index
        int atom1Idx = a_slices[sliceNum][atom1].getMatIdx();
        int atom1NumBasis = a_slices[sliceNum][atom1].getNumBasis();

        // Reset column offset
        offsetColumn = 0;

        // Iterate through secondary atoms
        for (int atom2=0; atom2<numAtoms2; atom2++){
          // Get atom 2 index
          int atom2Idx = a_slices[sliceNum+1][atom2].getMatIdx();
          int atom2NumBasis = a_slices[sliceNum+1][atom2].getNumBasis();

          // Loop over atom 1 basis functions
          for (int a1Basis=0; a1Basis<atom1NumBasis; a1Basis++){
            // Loop over atom 2 basis functions
            for (int a2Basis=0; a2Basis<atom2NumBasis; a2Basis++){

              /*
              std::cout << "TMPUMAT SIZE: " << arma::size(tmpUMat) << std::endl;
              std::cout << "HMAT SIZE: " << arma::size(a_HMat) << std::endl;
              std::cout << "TMPSUMAT SIZE: " << arma::size(tmpSUMat) << std::endl;
              std::cout << "SMAT SIZE: " << arma::size(a_SMat) << std::endl;

              std::cout << "UMAT ROW: " << offsetRow+a1Basis << std::endl;
              std::cout << "UMAT COLUMN: " << offsetColumn+a2Basis << std::endl;
              std::cout << "HMAT ROW: " << atom1Idx+a1Basis << std::endl;
              std::cout << "HMAT COLUMN: " << offsetColumn+a2Basis << std::endl;
              */

              // Set slice matrix elements
              tmpUMat(offsetRow+a1Basis,offsetColumn+a2Basis) = a_HMat(atom1Idx+a1Basis,atom2Idx+a2Basis);
              tmpSUMat(offsetRow+a1Basis,offsetColumn+a2Basis) = a_SMat(atom1Idx+a1Basis,atom2Idx+a2Basis);

              //std::cout << "SET MATRIX ELEMENTS" << std::endl;
            }
          }

          // Increment column
          offsetColumn += atom2NumBasis;
        }

        // Increment row
        offsetRow += atom1NumBasis;
      }

      // Add u matrices
      a_uSlices.push_back(tmpUMat);
      a_suSlices.push_back(tmpSUMat);
    }
  }
}
