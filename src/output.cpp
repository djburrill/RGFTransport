#include <output.h>

void output::writeData(std::string fileName, std::vector<double> a_c1, std::vector<double> a_c2){
  // Create stream file
  std::ofstream outFile;

  // Open file for writing
  outFile.open(fileName);

  // Write data
  for (int index=0; index<a_c1.size(); index++){
    outFile << a_c1[index] << " " << a_c2[index] << std::endl;
  }

  // Close file
  outFile.close();
}

void output::writeData(std::string fileName, std::vector<double> a_c1, std::vector<double> a_c2, std::vector<double> a_c3){
  // Create stream file
  std::ofstream outFile;

  // Open file for writing
  outFile.open(fileName);

  // Write data
  for (int index=0; index<a_c1.size(); index++){
    outFile << a_c1[index] << " " << a_c2[index] << " " << a_c3[index] << std::endl;
  }

  // Close file
  outFile.close();
}
