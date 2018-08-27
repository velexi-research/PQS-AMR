# setup-build.sh
#
# Script to set up build environment and directory.
#
# Usage:
#
# $ setup-build.sh

# --- Preparations

# Find top-level directory
TOP_DIR=`dirname "${BASH_SOURCE[0]}"`

# Set build directory
BUILD_DIR=$TOP_DIR/build

# Set SAMRAI directory
SAMRAI_DIR=$TOP_DIR/opt/SAMRAI

# --- Set compilers and flags

CC=gcc-8
CXX=g++-8
F77=gfortran-8
export CC CXX F77 

CXX_FLAGS="$CXX_FLAGS -std=c++11"

# --- Configure build

cmake . -B${BUILD_DIR} -DSAMRAI=${SAMRAI_DIR}

# --- Clean up

unset BUILD_DIR
unset SAMRAI_DIR
unset TOP_DIR
