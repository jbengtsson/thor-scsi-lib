
#lib_LIBRARIES          = libtracy.a
lib_LTLIBRARIES        = libtpsa.la libthor.la


## Todo
## conditional definition depending which library version is built
## or built both
#libtracy_a_SOURCES     = tracy_lib.cc
libtpsa_la_SOURCES     = tpsa_utils.cc field.cc
# libtpsa linear parts
libtpsa_la_SOURCES     += tpsa_lin.cc tpsa_lin_pm.cc
# libtpsa pm parts
# libtpsa_la_SOURCES     += tpsa_for_pm.cc # tpsa_for.cc
# libtest_la_SOURCES      = test_lib.cc
libthor_la_SOURCES     = radia2tracy.cc legacy/io.cc
libthor_la_SOURCES     += ety.cc eigenv.cc
libthor_la_SOURCES     += prtmfile.cc rdmfile.cc lsoc.cc
# libthor_la_SOURCES     += param.cc dynap.cc
libthor_la_SOURCES     +=  t2cell.cc t2lat.cc set_errors.cc
# Heavy template lifting below
# libthor_la_SOURCES     += t2elem.cc t2ring.cc thor_lib.cc
# libthor_la_SOURCES     += field.cc
libthor_la_SOURCES     += thor_lib.cc

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
AM_CXXFLAGS = -g -ggdb -O0  -pedantic -Wall  -std=gnu++14 -fPIC \
              -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
#AM_CXXFLAGS  += -pedantic -Wall  -Wcast-align -Wcast-qual -Wctor-dtor-privacy \
#                -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op \
#               -Wnoexcept -Wold-style-cast -Woverloaded-virtual \
#               -Wredundant-decls -Wshadow  -Wsign-promo \
#               -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef \
#               -Wno-unused
#AM_CXXFLAGS  += -Wmissing-include-dirs -Wsign-conversion -Wextra
AM_CXXFLAGS += -Werror
# not such a good choice fun with templates
AM_CXXFLAGS += -Wno-error=all
# better to catch mistakes with templates
AM_CXXFLAGS += -Wmissing-declarations
AM_CXXFLAGS += -Wno-error=unused-result -Wno-error=sign-compare -Wno-error=vla -Wno-error=old-style-cast -Wno-error=shadow -Wno-error=strict-overflow -Wno-error=sign-promo -Wno-error=format-nonliteral -Wno-error=switch-default
# AM_CXXFLAGS += -Wno-error=missing-include-dirs -Wno-error=old-style-cast -Wno-error=sign-conversion
# so that one fixes these
# AM_CXXFLAGS += -Werror
# of at least that one
AM_CXXFLAGS += -Werror=missing-declarations
# Fortran flags.
FFLAGS      = -g -O2 -Wall -fbounds-check -fdefault-integer-8 -mcmodel=large
