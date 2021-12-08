
#lib_LIBRARIES          = libtracy.a
lib_LTLIBRARIES        = libtpsa.la libthor.la


## Todo
## conditional definition depending which library version is built
## or built both
#libtracy_a_SOURCES     = tracy_lib.cc
libtpsa_la_SOURCES     = tpsa_utils.cc
libtpsa_la_SOURCES     += tpsa_lin.cc tpsa_lin_pm.cc
# libtpsa_la_SOURCES     += tpsa_for_pm.cc # tpsa_for.cc
libtest_la_SOURCES      = test_lib.cc
libthor_la_SOURCES     = radia2tracy.cc
libthor_la_SOURCES     +=  t2cell.cc t2lat.cc t2ring.cc t2elem.cc # thor_lib.cc

# AM_LDFLAGS flags:
#  {-static | -shared}
#  {-no-undefined | -undefined dynamic_lookup}
AM_LDFLAGS  = -shared

LIBS	    = -L../../TPSA \
              -lTPSALib -lLieLib \
              -lstdc++ \
              -lgsl \
              -L/usr/local/Cellar/armadillo/10.6.2/lib \
              -larmadillo
              $(GSL_LIBS)

# C Preprocessor flags.
AM_CPPFLAGS = -I../inc/  \
              -I/usr/local/Cellar/armadillo/10.6.2/include \
              -I/usr/local/Cellar/pybind11/2.7.1/include \
              -I/usr/local/Cellar/python@3.9/3.9.7/Frameworks/Python.framework/Versions/Current/include/python3.9 \
              $(GSL_CFLAGS)

# C++ flags:
#  {-g}                      debugging
#  {-pg}                     profiling
#  {-O0}                     optimization
#  {-Wall}                   warnings
#  {-Wno-error=all}          errors
#  {-fno-implicit-templates}
#  {-fPIC}
#  {-frepo}
AM_CXXFLAGS = -g -ggdb -O0 -Wall -Wno-error=all -std=gnu++14 -fPIC \
              -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC

# Fortran flags.
FFLAGS      = -g -O2 -Wall -fbounds-check -fdefault-integer-8 -mcmodel=large
