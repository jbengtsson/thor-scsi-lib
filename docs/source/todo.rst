Software update strategy
========================

1. Rework build process and compile files separately
2. Develop modern API for thor_scsi
3. Push code base to modern C++
4. Integrate thor_scsi as backbone for AT and PyAT
5. Implement python interface
6. Modernize reading lattice files
7. Use smart pointers for managing lattice elements
8. Make lattice managment transparant to python interface
9. Add solid set of unit tests
1. Review internal representation for cache friendlyness


Rework build process
--------------------

1. separate truncated power series in separate file
1. split up tracy.cc
2. rework remaining part in tracy.cc, which collects all 


