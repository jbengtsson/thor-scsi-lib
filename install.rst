Installing thor_scsi
====================


Dependencies
------------

- modern c++ compiler C++17 or better

    - std::shared_ptr
    - std::variant
    - std::ranges

- modern fortran compiler

- cmake
- pybind 11
- armadillo matrix library

- dependencies for building flame library:

   - flex and bison


Packages to be installed on Ubuntu 22 LTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ::

  sudo apt-get install bison flex cmake g++ gfortran libarmadillo-dev libboost-all-dev pybind11-dev python3-xarray



Checking out repository
-----------------------

.. ::

   git clone https://github.com/jbengtsson/thor-scsi-lib.git


change to the directory (persumably) `thor-scsi-lib`. If that was
successful please run

.. ::

   git submodule init
   git submodule update


Getting ready to build
----------------------

create a directory "build"

.. ::

   mkdir build


then change to this directory

.. ::

  cd build


then in this directory execute


.. ::

  cmake ..


This will create the build file. Typically this is a make file. In
case the cmake command fails, please remove at least the
`CMakeCache.txt` file in the build directory.

When cmake worked, trigger the build. In case you use make type

.. ::

  make

If build was successful use

.. ::

  cmake --install . --prefix=/path/to/install/to

with `/path/to/install/to` the absolute path of the directory you
would like to install to.

**NB: The libaries implementing the python interface will be
      currently installed in the source tree into directory
      `python/thor_scsi`. Have a look below for details
      of loading dynamic objects from non standard directories


Helping CMAKE find subcomponents
--------------------------------

Cmake checks that the version of required subcomponents is
sufficient. If it reports that one of the components is not
sufficiently new, I recommend to follow the following steps:

1. follow the instructions below required to make camke identify
   the component
2. After the cmake found the components  I recommend to

   1. remove the build directory
   2. create a new build directory
   3. run cmake in this clean directory.

Reason: cmake stores cache files and directories in the build
directory. These can still information from former cmake runs. In
my experience some rather strange configuration / build problems
are cured in this manner.



Up to date pybind11
~~~~~~~~~~~~~~~~~~~

If your version pybind 11 is rejected by cmake:

1. install it using pip

   ::
      pip3 install pybind11


   it can be that you have to use the `--user` flag so that it is
   installed within your environment.


2. help cmake find the installation. E.g. for a local installation
   on ubuntu (focal) it is typically found at

   ..highlight:: shell
      ls -d  $HOME/.local/lib/python3.8/site-packages/pybind11


   If still an too old version of pybind11 is found, please set
   the environment variable pybind11_DIR to the correct directory

   ..highlight:: shell
       export pybind11_DIR=$HOME/.local/lib/python3.8/site-packages/pybind11



Bison
-----

THe standard `bison` tool installed on mac os is not modern enough.
In our experience bison distributed with `brew` can be used. To
check if correct brew is installed in your shell run

::

    bison --config

The one installed on MAC OS is of major version 2 while version 3
is used for the parser used here. It seems that cmake does not
flag if the found bison binary is too old.

The following steps show what can be done, so that cmake will find
a sufficiently modern bison. So if not already installed, install
brew on your mac. Then follow `brew`  instruction to install
`bison`. Please find out where bison is located. (e.g.
`/usr/local/Cellar/bison/...`). Please add the directory of the
bison binary to the PATH variable (e.g. if you are using bash)


::

    export PATH=/path/to/bison:$PATH



Clear your build directory as explained above and check that a
sufficient modern bison version is found.

Loading dynamic objects from non standard locations
---------------------------------------------------

The libraries of thor-scsi-lib or the libraries for the python
interface can be installed in non standard places.

Linux
~~~~~
One solution can be to define the directory in LD_LIBRARY_PATH e.g.:

::
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/install/to/lib/





MAC OS
~~~~~~
One solution can be to define the directory in LD_LIBRARY_PATH e.g.:


::
    export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/install/to/lib/




Documentation
-------------

Requirements

* doxygen
* sphinx-doc
* breathe
* exhale
