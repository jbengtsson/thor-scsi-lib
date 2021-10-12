#!/bin/sh
CC=`pwd`/gcc_clang_stub.sh
export CC
# Will build but fail to link
/usr/local/Cellar/python\@3.9/3.9.6/bin/python3  setup.py build_ext -f

# link command
g++ -Wl,-no_compact_unwind -bundle -undefined dynamic_lookup -isysroot \
    /Library/Developer/CommandLineTools/SDKs/MacOSX10.14.sdk build/temp.macosx-10.14-x86_64-3.9/src/thor_py.o \
    -L../thor/src/.libs \
    -L../TPSA \
    -L/usr/local/Cellar/gsl/2.7/lib \
    -L/usr/local/lib \
    -L/usr/local/opt/openssl@1.1/lib \
    -L/usr/local/opt/sqlite/lib \
    -lthor \
    -lTPSALib -lLieLib \
    -lgsl  -lgslcblas \
    -o build/lib.macosx-10.14-x86_64-3.9/libthor.cpython-39-darwin.so \
    -mmacosx-version-min=10.14

# install it locally
/usr/local/Cellar/python\@3.9/3.9.6/bin/python3  setup.py build_ext -i
# if a n install is required
# /usr/local/Cellar/python\@3.9/3.9.6/bin/python3  setup.py install --user
