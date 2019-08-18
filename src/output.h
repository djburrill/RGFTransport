/** Output/Print Methods.
*/
#ifndef _output_h_
#define _output_h_
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

class output{
public:
  /** Default constructor
  */
  output(){};

  /** Write two column data
  */
  void writeData(std::string fileName, std::vector<double> a_c1, std::vector<double> a_c2);

  /** Write three column data
  */
  void writeData(std::string fileName, std::vector<double> a_c1, std::vector<double> a_c2, std::vector<double> a_c3);

private:

};
#endif
