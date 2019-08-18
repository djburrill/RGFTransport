/** Atom container.
*/
#ifndef _atom_h_
#define _atom_h_
#include <vector>
#include <iostream>

class atom{
public:
  /** Set up atom
  */
  atom(std::string a_label, double a_x, double a_y, double a_z);

  /** Get atomic position
  */
  std::vector<double> getPos();

  /** Get atomic label
  */
  std::string getLabel();

  /** Print atomic information
  */
  void print();

  /** Set atom number
  */
  void setAtomNum(int a_atomNum);

  /** Get atom number
  */
  int getAtomNum();

  /** Set number of basis functions
  */
  void setNumBasis(int a_basisNum);

  /** Get number of basis function
  */
  int getNumBasis();

  /** Set matrix index
  */
  void setMatIdx(int a_matIdx);

  /** Get matrix index
  */
  int getMatIdx();

private:
  // Private variables
  std::vector<double> m_position;
  std::string m_label;
  int m_atomNum;
  int m_numBasis;
  int m_matIdx;                     /**< First index of atom basis in Hamiltonian matrix */
};
#endif
