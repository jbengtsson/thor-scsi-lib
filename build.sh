#!/bin/sh

# Script to build thor-scsi-lib.

home_dir=`pwd`
echo "home_dir =" $home_dir

# Clone repository & submodules - clone leaves submodules empty.
if false; then
    git clone --recursive https://github.com/jbengtsson/thor-scsi-lib.git
    git submodule update --init --recursive
fi

# Pull all changes in the repository:
if false; then
    git pull --recurse-submodules
fi

# Checkout local branches for gtpsa from CERN.
if ! false; then
    cd src/gtpsa
    git checkout gtpsa_jb
    cd mad-ng
    git checkout 11fcbfeb
fi

# Create the build directory.
if ! false; then
    mkdir build
fi

# Make & install the thor-scsi & gtpsa libraries.
if ! false; then
    cd $home_dir/build
    cmake ..
    make -j8
    cmake --install . --prefix ../local
fi

# Validate the thor-scsi & gtpsa libraries.
if false; then
    cd $home_dir/build
    make test
fi

# Build the Python interfaces.
if ! false; then
    export THOR_SCSI_LIB=$home_dir
    export gtpsa_PREFIX=$THOR_SCSI_LIB/local
    export thor_scsi_PREFIX=$THOR_SCSI_LIB/local
    echo "\n\$THOR_SCSI_LIB set to: " $THOR_SCSI_LIB
    echo "\$gtpsa_PREFIX set to:      " $gtpsa_PREFIX
    echo "\$thor_scsi_PREFIX set to:  " $thor_scsi_PREFIX
    cd $home_dir/src/gtpsa/python
    pip install .
    cd $home_dir/python
    pip install .
fi

# Create a local Python environment.
if false; then
    python -m venv $home_dir
    source $home_dir/bin/activate
fi
