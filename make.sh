#!/bin/sh

# Script to build thor-scsi-lib.
# Cd to local directory for repository before running script.

home_dir=`pwd`
echo "home_dir =" $home_dir

# Clone repository & submodules - clone by default leaves submodules empty.
if false; then
    git clone --recursive https://github.com/jbengtsson/thor-scsi-lib.git
    # git submodule update --init --recursive
    # git pull --recurse-submodules
fi

# Create the build directory.
if false; then
    cd $home_dir
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

# Create a local Python environment.
if false; then
    python3 -m venv $home_dir/venv
    # Source vs. "." source only works for bash and not sh.
    . $home_dir/venv/bin/activate
    # To terminate: deactivate.
fi

# Install required libraries.
if  false; then
    # Upgrade from pip-24.0 to pip-24.1.
    pip install --upgrade pip
    # Upgrade: wheel, setuptools, and pip.
    pip install wheel setuptools pip --upgrade
    pip install pybind11

    pip install numpy
    pip install scipy
    pip install xarray
    pip install matplotlib
fi

# Build the Python interfaces.
if ! false; then
    export THOR_SCSI_LIB=$home_dir
    export gtpsa_PREFIX=$THOR_SCSI_LIB/local
    export thor_scsi_PREFIX=$THOR_SCSI_LIB/local
    echo "\n\$THOR_SCSI_LIB set to:     " $THOR_SCSI_LIB
    echo "\$gtpsa_PREFIX set to:      " $gtpsa_PREFIX
    echo "\$thor_scsi_PREFIX set to:  " $thor_scsi_PREFIX
    cd $home_dir/src/gtpsa/python
    pip install .
    cd $home_dir/python
    pip install .
fi

# Export path for libraries.
if false; then
    # Linux.
    echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$home_dir/local/lib" \
	 >> ~/.bashrc
    # Macbook.
    # echo "export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$home_dir/local/lib" \
    # 	 >> ~/.bashrc
fi
