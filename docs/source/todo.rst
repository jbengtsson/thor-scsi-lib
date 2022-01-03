Software update strategy
========================

1. Rework build process and compile files separately
2. Develop modern API for thor_scsi
3. Push code base to modern C++
  a. namespaces
  b. extensive use of standard library
  c. Library for arrays matrix and vectors: use library for parallel execution
4. Integrate thor_scsi as backbone for AT and PyAT
5. Implement python interface based on pybind11
6. Modernize reading lattice files
7. Use smart pointers for managing lattice elements
8. Make lattice managment transparant to python interface
9. Add solid set of unit tests
10. Review internal representation for cache friendlyness
11. Calculations (e.g. integration steps) separate from memory management


Rework build process
--------------------

1. separate truncated power series in separate file
2. split up tracy.cc
3. rework remaining part in tracy.cc, which collects all 


