# RGFTransport

## Description
The non-equilibrium Green's function (NEGF) method is widely used among electronic structure codes to calculate electron transport through a device bridging a number of 

## Dependencies
- Armadillo Linear Algebra Library [http://arma.sourceforge.net]

## Build
- Navigate to Code/exec
- Run 'make all'

### Special Considerations
- I do not think there is anything that needs to be changed in the make file. Just make sure the armadillo library is able to be found using your PATH environment variable.
- Use the utils as a guide for how to perform analysis. Some of them are written for debugging purposes so the expected file formatting may be different from what the program outputs. Check to be sure!
