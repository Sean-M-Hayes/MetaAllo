# MetaAllo
MetaAllo is a platform for the simulation and analysis of an Allometric Metacommunity model. 

The code here additionally requires an ODE solver to run. As written it expects to use the FORTRAN77 ODEPACK solver LSODA (https://computation.llnl.gov/casc/odepack/), although care must be taken to make sure that FORTRAN Common blocks do not cause issues with the multithreading present in the C++ code. 

A project file (Metacommunity_Model_UI.pro) used to compile the program with Qt Creator is included and is useful for determining requirements for compiling. In addition, the program depends on many Boost libraries, the OpenMP API for multithreading, and the Eigen header library. 

This code was written and revised sporadically as a first project for learning C++, and is tragically naive to commenting ettiquite and formatting standards. 
