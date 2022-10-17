# Build command on MacBook:
#
#   unset CC
#   unset CXX
#   CC=`pwd`/gcc_clang_stub.sh python3 setup.py build_ext -i

import logging
import os
import sys
from setuptools import setup
from distutils.file_util import copy_file

from pybind11.setup_helpers import Pybind11Extension, build_ext

import gsl_conf
import sys

# Command line:
#   \rm -rf build
#   CC=g++-11 python3 setup.py <build|build_ext|install|--help>

# Make a copy of README.rst ... required by setup.cfg
t_dir = os.path.dirname(__file__)
par_dir = os.path.normpath(os.path.join(t_dir, os.pardir))
readme_name = "README.rst"
copy_file(
    os.path.join(par_dir, readme_name), os.path.join(t_dir, readme_name), update=True
)

d = gsl_conf.gsl_config()

# How to define where the thor scsi library is located?
# here there are some examples
prefix = os.path.abspath(os.path.join(os.path.dirname(__name__), os.pardir, os.pardir))
prefix = os.path.abspath(os.path.join(os.environ["HOME"], ".local"))
prefix = os.path.abspath(os.path.join(os.path.dirname(__name__), os.pardir, "local"))
# prefix = os.path.abspath(os.path.join(os.environ["HOME"], ".local"))

boost_prefix="/usr/include"
if sys.platform == "darwin":
    boost_prefix=os.path.join("/", "usr", "local", "include")

from pybind11.setup_helpers import ParallelCompile

# Optional multithreaded build
try:
    os.environ["NPY_NUM_BUILD_JOBS"]
except KeyError:
    os.environ["NPY_NUM_BUILD_JOBS"] = "5"

ParallelCompile("NPY_NUM_BUILD_JOBS").install()

ext_modules = [
    # Pybind11Extension(
    #    "testbind",
    #    sorted(['src/test_bind.cc', 'src/mytest.cc']),
    #    include_dirs=['.'],
    # ),
    Pybind11Extension(
        "flame",
        ["src/flame.cc"],
        include_dirs=[d["gsl_include"]] + [os.path.join(prefix, "include"), boost_prefix],
        library_dirs=([os.path.join(prefix, "lib")]),
        libraries=["flame", "flame_core"],
    ),
    Pybind11Extension(
        "lib",
        sorted(
            [
                "src/thor_scsi.cc",
                "src/arma.cc",
                "src/tps.cc",
                "src/config_type.cc",
                "src/elements.cc",
                "src/accelerator.cc",
            ]
        ),
        include_dirs=[d["gsl_include"]] + [os.path.join(prefix, "include"), boost_prefix],
        # define_macros=[("_GLIBCXX_DEBUG", 1), ("_GLIBCXX_DEBUG_PEDANTIC", 1)],
        library_dirs=(
            [os.path.join(prefix, "lib")]
            # ["../../engine/lib"]
            + [d["gsl_lib_dir"]]
        ),
        libraries=["thor_scsi", "thor_scsi_core", "tpsa_lin", "flame", "flame_core"]
        + d["gsl_libs"],
    ),
]


setup(
    cmdclass={"build_ext": build_ext}, ext_package="thor_scsi", ext_modules=ext_modules
)
