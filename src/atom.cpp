#include <atom.h>

atom::atom(std::string a_label, double a_x, double a_y, double a_z){
  // Set label and position
  m_label = a_label;
  m_position.push_back(a_x);
  m_position.push_back(a_y);
  m_position.push_back(a_z);
}

std::vector<double> atom::getPos(){
  return m_position;
}

std::string atom::getLabel(){
  return m_label;
}

void atom::print(){
  std::cout << m_label << ' ' << m_position[0] << ' ' << m_position[1] << ' ' << m_position[2] << std::endl;
}

void atom::setAtomNum(int a_atomNum){
  m_atomNum = a_atomNum;
}

int atom::getAtomNum(){
  return m_atomNum;
}

void atom::setNumBasis(int a_numBasis){
  m_numBasis = a_numBasis;
}

int atom::getNumBasis(){
  return m_numBasis;
}

void atom::setMatIdx(int a_matIdx){
  m_matIdx = a_matIdx;
}

int atom::getMatIdx(){
  return m_matIdx;
}
