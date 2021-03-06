Introduction
-----------------------------------------------
This guide is a developers guide to the CCPi Quantification project (If you are looking for how to use the CCPi Quantification plugins then please read "CCPi Quantification user guide"). CCPi Quantification has set of algorithms that are re-engineered from the code/algorithms developed by the X-Ray imaging community. 
The aim of the project is to make these scientific algorithms/code easily accessible, runs efficiently and better maintainable.  

* accessibility: is achieved through making the algorithms available as plugins in popular user applications. 
* efficiency: this is achieved through reengineering the code to run in parallel and using faster mathematical library.
* maintainability: the code is available in git repository, includes doxygen comments and integrated with continious integration for testing.
   
This project currently supports plugin for three popular applications.

* Avizo    (Versions: 7.1.1, 8.x, 9.0, 9.1.1, 9.2)
* Paraview (Versions: 5.2)
* ImageJ   (Versions: 150)
