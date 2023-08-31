#!/bin/bash

set -ex

# Compile options

# Default (optimized) settings.  Adept's fast exponential gives
# slightly different results - actually better in terms of the
# smoothness of the performance as a function of number of streams -
# but for repeatability we prefer the standard exponential
#CXXFLAGS="-Wall -g -O2 -march=native -std=c++11 -DADEPT_FAST_EXPONENTIAL"

# The values in the paper do not use the fast exponential
CXXFLAGS="-Wall -g -O2 -march=native -std=c++11"

# Optimized but with checking
#CXXFLAGS="-Wall -g -O2 -march=native -std=c++11 -DADEPT_BOUNDS_CHECKING -DADEPT_INIT_REAL_SNAN"

# Debug settings
#CXXFLAGS="-Wall -g -O0 -march=native -std=c++11 -DADEPT_BOUNDS_CHECKING -DADEPT_INIT_REAL_SNAN"

# Location of Adept automatic differentiation library
ADEPT_VER=adept-2.1.2
ADEPT_DIR=/home/parr/apps/$ADEPT_VER
ADEPT_FLAGS="--with-adept=$ADEPT_DIR"

# Location of NetCDF-4 library
module load netcdf4
NETCDF_FLAGS="--with-netcdf=$NETCDF4_DIR"

# Set install location
INSTALL_DIR=/var/tmp/$USER/optimize_angles

# Call configure script
./configure --prefix "$INSTALL_DIR" "CXXFLAGS=$CXXFLAGS" $ADEPT_FLAGS $NETCDF_FLAGS "LDFLAGS=$LDFLAGS" $@
