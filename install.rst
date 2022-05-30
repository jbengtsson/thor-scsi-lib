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

   git clone


change to the directory then


for the time being run:

..::

   git checkout flame-integration


if that was successful please run

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


This willl create the build file. Typicslly this is a make file. In case you use make
type

.. ::

  make

If build was successful use

.. ::

  cmake --install . --prefix=/path/to/install/to

with `/path/to/install/to` the absolute path of the directory you would like to install to.


Helping CMAKE find subcomponents
--------------------------------

Cmake checks that the version of required subcomponents is sufficient. If it reports that one of the components
is not sufficiently new, I recommend to follow the following steps:

1. follow the instructions below required to make camke identify the component
2. After the cmake found the components  I recommend to

   1. remove the build directory
   2. create a new build directory
   3. run cmake in this clean directory.

Reason: cmake stores cache files and directories in the build directory. These can still information from former cmake
runs. In my experience some rather strange configuration / build problems are cured in this manner



Up to date pybind11
~~~~~~~~~~~~~~~~~~~

If your version pybind 11 is rejected by cmake:

1. install it using pip

   ```python

   pip3 install pybind11
   ```

   it can be that you have to use the `--user` flag so that it is installed
   within your environment.


2. help cmake find the installation. E.g. for a local installation on ubuntu (focal)
   it is typically found at

   ```shell

    ls -d  $HOME/.local/lib/python3.8/site-packages/pybind11

   ```

   If still an too old version of pybind11 is found, please set the environment
   variable pybind11_DIR to the correct directory
   ```shell
    export pybind11_DIR=$HOME/.local/lib/python3.8/site-packages/pybind11
   ```


Bison
-----

THe standard `bison` tool installed on mac os is not modern enough In our experience bison distributed
with `brew` can be used. To check if correct brew is installed in your shell run

.. ::

   bison --config

THe one installed on MAC OS is of major version 2 while version 3 is used for the parser used here.


So if not already installed, install brew on your mac. Then follow `brew`  instruction to install `bison`. Now the
`PATH` variable needs to be modifed so that cmake will find bison. Now you need to find where `bison`




Documentation
-------------

Requirements

* doxygen
* sphinx-doc
* breathe
* exhale
