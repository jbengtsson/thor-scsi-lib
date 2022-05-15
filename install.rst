Installing thor_scsi
====================

Dependencies
------------

- modern c++ compiler C++17 or better

    - std::shared_ptr
    - std::variant
    - std::ranges

- cmake
- pybind 11


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
