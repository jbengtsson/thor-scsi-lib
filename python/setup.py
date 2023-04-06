# Build command on MacBook:
#
#   unset CC
#   unset CXX
#   CC=`pwd`/gcc_clang_stub.sh python3 setup.py build_ext -i

import logging
import os
import sys
from setuptools import setup, find_packages
from distutils.file_util import copy_file

from pybind11.setup_helpers import Pybind11Extension, build_ext

import gsl_conf
import sys

logger = logging.getLogger("thor-scsi-build")
# Command line:
#   \rm -rf build
#   CC=g++-11 python3 setup.py <build|build_ext|install|--help>

# Make a copy of README.rst ... required by setup.cfg
# only valid if setup.py is called within the
# thor scsi tree
t_dir = os.path.dirname(__file__)
par_dir = os.path.normpath(os.path.join(t_dir, os.pardir))
readme_name = "README.rst"
readme_path = os.path.join(par_dir, readme_name)

do_readme_copy = True
try:
    os.stat(readme_path)
except FileNotFoundError:
    do_readme_copy = False

if do_readme_copy:
    copy_file(readme_path, os.path.join(t_dir, readme_name), update=True)

d = gsl_conf.gsl_config()

prefix = None

# How to define where the thor scsi library is located?
# Some hard coded preferences for your convienience
# prefix = os.path.abspath(os.path.join(os.path.dirname(__name__), os.pardir, os.pardir))
# prefix = os.path.abspath(os.path.join(os.environ["HOME"], ".local"))

# deliberatley chosen _PREFIX instead of DIR ...
# not to interfear with cmake's variables
try:
    prefix = os.environ["thor_scsi_PREFIX"]
except KeyError as ke:
    logger.info(f"no environment variable thor_scsi_PREFIX: ke")

# a hack for mac
boost_prefix = "/usr/include"
if sys.platform == "darwin":
    boost_prefix = os.path.join("/", "usr", "local", "include")

from pybind11.setup_helpers import ParallelCompile

# Optional multithreaded build
try:
    os.environ["NPY_NUM_BUILD_JOBS"]
except KeyError:
    os.environ["NPY_NUM_BUILD_JOBS"] = "5"

ParallelCompile("NPY_NUM_BUILD_JOBS").install()

# Build up include and lib directories respecting prefix
include_dirs = [boost_prefix]
library_dirs = []
if prefix:
    include_dirs += [os.path.join(prefix, "include")]
    library_dirs += [os.path.join(prefix, "lib")]

ext_modules = [
    # Pybind11Extension(
    #    "testbind",
    #    sorted(['src/test_bind.cc', 'src/mytest.cc']),
    #    include_dirs=['.'],
    # ),
    Pybind11Extension(
        "pyflame",
        ["src/flame.cc"],
        include_dirs=include_dirs,
        libraries=["flame", "flame_core"],
        library_dirs=library_dirs,
    ),
    Pybind11Extension(
        "lib",
        sorted(
            [
                "src/accelerator.cc",
                "src/aperture.cc",
                "src/config_type.cc",
                "src/elements.cc",
                "src/enums.cc",
                "src/flame.cc",
                "src/interpolation.cc",
                "src/observer.cc",
                "src/pybind_test.cc",
                "src/radiation.cc",
                "src/tps.cc",
                "src/thor_scsi.cc",
                "src/custom.cc"
            ]
        ),
        # Required for MacBook llvm C++ compiler.
        define_macros=[("GTPSA_DEFINE_BOOL", 1)],
        include_dirs=[d["gsl_include"]] + include_dirs,
        # define_macros=[("_GLIBCXX_DEBUG", 1), ("_GLIBCXX_DEBUG_PEDANTIC", 1)],
        library_dirs=(
            library_dirs
            # ["../../engine/lib"]
            + [d["gsl_lib_dir"]]
        ),
        libraries=[
            "thor_scsi_gtpsa",
            "thor_scsi_core_gtpsa",
            "gtpsa",
            "tpsa_lin",
            "flame",
            "flame_core",
        ]
        + d["gsl_libs"],
    ),
]


setup(
    cmdclass={"build_ext": build_ext},
    ext_package="thor_scsi",
    ext_modules=ext_modules,
    packages=find_packages(),
)
