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


change to the directory then run

.. ::

   git submodule init
   git submodule update



Helping CMAKE find subcomponents
--------------------------------

Cmake checks that the version of required subcomponents is sufficient.

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

Documentation
-------------

Requirements

* doxygen
* sphinx-doc
* breathe
* exhale
