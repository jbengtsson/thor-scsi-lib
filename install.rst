.. _install.rst:

Installing thor-scsi-lib
========================

`thor-scsi-lib` consists of a c++ library and python wrapper to it. Here
it is first described how to build the library and then how to build the
python wrapper.

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

- modern python3

Packages to be installed on Ubuntu 22 LTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install packages that could be missing for the c++ library

.. code:: shell

  sudo apt-get install bison flex cmake g++ gfortran libarmadillo-dev libboost-all-dev


The following packages could be missing for the python wrapper

.. code:: shell

   sudo apt-get install bison flex cmake g++ gfortran libarmadillo-dev libboost-all-dev pybind11-dev python3-xarray pybind11-dev python3-xarray


Checking out repository
-----------------------

First clone the repository using

.. code:: shell

   git clone https://github.com/jbengtsson/thor-scsi-lib.git


change to the directory (persumably) `thor-scsi-lib`.

Then initialise submodules using the following command

.. code:: shell

   git submodule update --init --recursive

*NB*: this command currently will pull a subrepository (`cmake4epics`).
This repository currently does not support (llvm/clang). Thus build on
MAC currently fails. A fix is currently worked on.



Getting ready to build
----------------------

create a directory "build"

.. code:: shell

   mkdir build


then change to this directory

.. code:: shell

  cd build


then in this directory execute

.. code:: shell

  cmake ..


This will create the build file. Typically this is a make file. In
case the cmake command fails, please remove at least the
`CMakeCache.txt` file in the build directory. If this steps fails,
find some hints how to solve them in section
"Helping CMAKE find subcomponents" :ref:`cmake-find-subcomponents`.


When cmake worked, trigger the build. In case you use `make` type

.. code:: shell

  make


The build can be verified executing the tests using

.. code:: shell

   make test


If build was successful use

.. code:: shell

  cmake --install . --prefix < path to install to e.g.: > ../local


with `path to install to` the absolute path of the directory you
would like to install to.

**NB**: The libaries implementing the python interface will be
currently installed in the source tree into directory
`python/thor_scsi` and src/gtpsa/python.
Have a look below for details
of loading dynamic objects from non standard directories
if you want to use these. The python wrapper and module
can be installed using `setup.py` too.



Installing python module thor_scsi and gtpsa
--------------------------------------------

Currently the python wrapper is automatically built when the c++ library is built.
Additionally a `setup.py` script is provided that can be used to use the standard
python install procedure.

Before you can use this script, you need to build the c++ library and install it
 to some path (called `/path/to/install/to` above).

Directories with python modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two python modules are provided

* gtpsa: directory src/gtpsa/python
* thor_scsi: directory python/

Recommandation is to first build gtpsa and then thor scsi.
The description below refers to both of them. Both directories are 
refered to as `python` directory below.

Installation instruction for one of the packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The process of building the python package depends on the C++ headers and libraries.
Thus it needs to know where these are installed. The process described below will
use the directory provided by the corresponding environment variables:

* gtpsa_PREFIX for the `gtpsa` package
* thor_scsi_PREFIX for the `thor_scsi` package.

Please note that when `thor_scsi` is built `pyflame` is built too.

Change into the repositories `python` directory. Edit the
`setup.py` file and define the variable `prefix` to contain the path you installed
the C++ library to.

For `gtpsa` this would be

.. code:: shell

    export gtpsa_PREFIX=path/to/install/to

For `thor_scsi` this would be

.. code:: shell

    export thor_scsi_PREFIX=path/to/install/to


As soon that has been done, you should be able to use e.g.

.. code:: shell

   python setup.py build


to build the module and

.. code:: shell

   python setup.py install


to install the module.


Alternatively you could use `pip` e.g.

.. code:: shell

   pip install .

to install the package.

If you are curios how and where to get pip see the link below.

https://pip.pypa.io/en/stable/installation/

.. _cmake-find-subcomponents:

Helping CMAKE find subcomponents
--------------------------------

Here some information if cmake above fails.

Cmake checks that the version of required subcomponents is
sufficient. If it reports that one of the components is not
sufficiently new, I recommend to follow the following steps:

1. follow the instructions below required to make cmake identify
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

   .. code:: shell

      pip3 install pybind11


   it can be that you have to use the `--user` flag so that it is
   installed within your environment.


2. help cmake find the installation. E.g. for a local installation
   on ubuntu (focal) it is typically found at

   .. code:: shell

      ls -d  $HOME/.local/lib/python3.8/site-packages/pybind11


   If still an too old version of pybind11 is found, please set
   the environment variable pybind11_DIR to the correct directory
      e.g.

   .. code:: shell

       export pybind11_DIR=$HOME/.local/lib/python3.8/site-packages/pybind11



Bison
-----

THe standard `bison` tool installed on mac os is not modern enough.
In our experience bison distributed with `brew` can be used. To
check if correct brew is installed in your shell run

.. code:: shell

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


.. code:: shell

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

.. code:: shell

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/install/to/lib/


macOS
~~~~~~
One solution can be to define the directory in LD_LIBRARY_PATH e.g.:


.. code:: shell

    export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/path/to/install/to/lib/


Python 3.12 is required to install SciPy via wheels.
Python 3.13+ currently requires building SciPy from source and is
not supported due to OpenMP toolchain limitations.


Quick sanity check
~~~~~~~~~~~~~~~~~~

python - <<'EOF'
import numpy, scipy, matplotlib, xarray
print("NumPy:", numpy.__version__)
print("SciPy:", scipy.__version__)
EOF


Documentation
-------------

Requirements

* doxygen
* sphinx-doc
* breathe
* exhale
