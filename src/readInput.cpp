#include <readInput.h>

readInput::readInput(std::string a_inFileName){
  // Variables
  int idxCounter = 0;
  int lOffset = 0;

  // Read input file
  readInputFile(a_inFileName);

  // Set left lead offset
  //lOffset = m_lLead.size();
  lOffset = 0;

  // Set atom parameters
  // Left lead
  for (int i=0; i<m_lLead.size(); i++){
    // Atom number
    m_lLead[i].setAtomNum(i);

    // Number of basis functions
    m_lLead[i].setNumBasis(m_numBasis[m_lLead[i].getLabel()]);
  }

  // Device
  for (int i=0; i<m_molecule.size(); i++){
    // Atom number
    m_molecule[i].setAtomNum(lOffset+i);

    // Number of basis functions
    m_molecule[i].setNumBasis(m_numBasis[m_molecule[i].getLabel()]);

    // Matrix index
    m_molecule[i].setMatIdx(lOffset+idxCounter);
    idxCounter += m_molecule[i].getNumBasis();
  }

  // Right lead
  for (int i=0; i<m_rLead.size(); i++){
    // Atom number
    m_rLead[i].setAtomNum(i);

    // Number of basis functions
    m_rLead[i].setNumBasis(m_numBasis[m_rLead[i].getLabel()]);
  }

  // Calculate the number of basis functions per unit cell in the leads
  calcBasisUCLead();

  // Print basis
  print_Basis();
}

void readInput::calcBasisUCLead(){
  // Left lead
  for (int i=0; i<m_llucAtoms; i++){
    m_numBasisUCL += m_numBasis[m_lLead[i].getLabel()];
  }

  // Right lead
  for (int i=0; i<m_rlucAtoms; i++){
    m_numBasisUCR += m_numBasis[m_rLead[i].getLabel()];
  }
}

void readInput::print_Basis(){
  for (auto const &entry : m_numBasis){
    auto const &key = entry.first;
    //std::cout << key << " : " << m_numBasis[key] << std::endl;
  }
}

void readInput::readInputFile(std::string a_inFileName){
  // Variables
  std::ifstream inFile(a_inFileName);
  std::string fileLine;
  std::string tmpString;

  // Open file
  if (inFile.is_open()){
    // Parse file
    while (getline(inFile,fileLine)){
      // Trim line
      tmpString = trimString(fileLine);

      // Skip blank lines
      if (tmpString.size() == 0){
        continue;
      }

      // Check block name
      if (tmpString[0] == '/'){
        m_inputCard = tmpString.substr(1);
        continue;
      }

      // Read properties
      if (m_inputCard == "prop"){
        readProperty(tmpString);
      }
      // Read device atoms
      else if (m_inputCard == "device"){
        readAtomPos(m_molecule,tmpString);
      }
      // Read left lead atoms
      else if (m_inputCard == "llead"){
        readAtomPos(m_lLead,tmpString);
      }
      // Read right lead atoms
      else if (m_inputCard == "rlead"){
        readAtomPos(m_rLead,tmpString);
      }
      // Read right lead atoms
      else if (m_inputCard == "basis"){
        readBasis(tmpString);
      }
      // Unknown block
      else{
        std::cout << "ERROR IN INPUT FILE. UNKNOWN BLOCK." << std::endl;
        return;
      }
    }

    // Close file
    inFile.close();
  }
  // Return error message if cannot open file
  else{
    std::cout << "Unable to open input file: " << a_inFileName << std::endl;;
  }
}

void readInput::readProperty(std::string a_inString){
  // Variables
  std::size_t splitPos;
  std::string tag;
  std::string value;

  // Find position of '='
  splitPos = a_inString.find('=');

  // Split line into tag and value
  tag = a_inString.substr(0,splitPos-1);
  value = a_inString.substr(splitPos+1);

  // Trim strings
  tag = trimString(tag);
  value = trimString(value);

  // Set property value
  if (tag == "minSliceDist"){
    m_minSliceDist = std::stod(value);
  }
  else if (tag == "llucatoms"){
    m_llucAtoms = std::stoi(value);
  }
  else if (tag == "rlucatoms"){
    m_rlucAtoms = std::stoi(value);
  }
  else if (tag == "lhfilename"){
    readMatrix(m_lhMat,value);
  }
  else if (tag == "lsfilename"){
    readMatrix(m_lsMat,value);
  }
  else if (tag == "rhfilename"){
    readMatrix(m_rhMat,value);
  }
  else if (tag == "rsfilename"){
    readMatrix(m_rsMat,value);
  }
  else if (tag == "dhfilename"){
    readMatrix(m_dhMat,value);
  }
  else if (tag == "dsfilename"){
    readMatrix(m_dsMat,value);
  }
  else if (tag == "imag_eta"){
    m_imag_eta = std::stod(value);
  }
  else if (tag == "sc_conv"){
    m_sc_conv = std::stod(value);
  }
  else if (tag == "vlow"){
    m_vlow = std::stod(value);
  }
  else if (tag == "vhigh"){
    m_vhigh = std::stod(value);
  }
  else if (tag == "dv"){
    m_dv = std::stod(value);
  }
  else if (tag == "elow"){
    m_elow = std::stod(value);
  }
  else if (tag == "ehigh"){
    m_ehigh = std::stod(value);
  }
  else if (tag == "de"){
    m_de = std::stod(value);
  }
  else{
    std::cout << "UNKNOWN TAG: " << tag << std::endl;
    exit(EXIT_FAILURE);
  }
}

void readInput::readAtomPos(std::vector<atom>& a_molecule, std::string a_inString){
  // Variables
  std::size_t splitPos;
  std::string tmpString;
  std::string labelStr;
  std::string xStr;
  std::string yStr;
  std::string zStr;

  // Trim string
  tmpString = trimString(a_inString);

  // Split at first space
  splitPos = tmpString.find(' ');
  labelStr = tmpString.substr(0,splitPos);
  labelStr = trimString(labelStr);
  tmpString = trimString(tmpString.substr(splitPos));

  // Split at second space
  splitPos = tmpString.find(' ');
  xStr = tmpString.substr(0,splitPos);
  xStr = trimString(xStr);
  tmpString = trimString(tmpString.substr(splitPos));

  // Split at third space
  splitPos = tmpString.find(' ');
  yStr = tmpString.substr(0,splitPos);
  yStr = trimString(yStr);
  tmpString = trimString(tmpString.substr(splitPos));

  // Set Z string
  zStr = tmpString;

  // Assign atom properties
  atom tmpAtom(labelStr,std::stod(xStr),std::stod(yStr),std::stod(zStr));
  a_molecule.push_back(tmpAtom);
}

void readInput::readBasis(std::string a_inString){
  // Variables
  std::size_t splitPos;
  std::string tmpString;
  std::string label;
  std::string value;

  // Trim string
  tmpString = trimString(a_inString);

  // Find position of ' '
  splitPos = a_inString.find(' ');

  // Split line into label and value
  label = a_inString.substr(0,splitPos);
  value = a_inString.substr(splitPos);

  // Trim label and value strings
  label = trimString(label);
  value = trimString(value);

  // Map label and value
  m_numBasis[label] = std::stoi(value);
}

void readInput::readMatrix(arma::mat& a_inMat, std::string a_fileName){
  a_inMat.load(a_fileName);
}

std::string readInput::trimStringLeft(std::string inString, char trimChar){
  // Variables
  int counter = 0;

  // Iterate through string to remove trimChar
  while (counter < inString.size()){
    if (inString[counter] == trimChar){
      counter += 1;
    }
    else{
      break;
    }
  }

  // Return conditions
  if (counter == inString.size()){
    return "";
  }
  else{
    return inString.substr(counter);
  }
}

std::string readInput::trimStringLeft(std::string inString){
  // Variables
  int counter = 0;

  // Iterate through string to remove trimChar
  while (counter < inString.size()){
    if (inString[counter] == ' '){
      counter += 1;
    }
    else{
      break;
    }
  }

  // Return conditions
  if (counter == inString.size()){
    return "";
  }
  else{
    return inString.substr(counter);
  }
}

std::string readInput::trimStringRight(std::string inString, char trimChar){
  // Variables
  int counter = inString.size()-1;
  std::string tmpString = inString;

  // Iterate through string to remove trimChar
  while (counter > -1){
    if (inString[counter] == trimChar){
      counter -= 1;
    }
    else{
      break;
    }
  }

  // Return conditions
  if (counter == -1 ){
    return "";
  }
  else{
    tmpString.resize(counter+1);
    return tmpString;
  }
}

std::string readInput::trimStringRight(std::string inString){
  // Variables
  int counter = inString.size()-1;
  std::string tmpString = inString;

  // Iterate through string to remove trimChar
  while (counter > -1){
    if (inString[counter] == ' '){
      counter -= 1;
    }
    else{
      break;
    }
  }

  // Return conditions
  if (counter == -1 ){
    return "";
  }
  else{
    tmpString.resize(counter+1);
    return tmpString;
  }
}

std::string readInput::trimString(std::string inString, char trimChar){
  // Variables
  std::string tmpString = inString;

  // Trim left
  tmpString = trimStringLeft(tmpString, trimChar);

  // Trim right
  tmpString = trimStringRight(tmpString, trimChar);

  return tmpString;
}

std::string readInput::trimString(std::string inString){
  // Variables
  std::string tmpString = inString;

  // Trim left
  tmpString = trimStringLeft(tmpString);

  // Trim right
  tmpString = trimStringRight(tmpString);

  return tmpString;
}
