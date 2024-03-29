#!/bin/sh
#
# Specify (M1, M2, mu, tanb) 
# candidate points.
# 
# Compile and link the code that decomposes and diagonalizes 
# the neutralino and chargino matrices.
# This is underlying Fortran wrapped in C++.
# It uses LAPACK (and BLAS) and ROOT
#
echo 'Check ROOT environment variables are set up OK'
echo ' '
echo $ROOTSYS
echo ' '
echo $PATH
echo ' '
echo $LD_LIBRARY_PATH

# Since we are terribly inefficient turn off printing details !
# (do this using ewdriver ...)
ln -s lprint_true.inc lprint.inc

gfortran -c electroweakino.f90

g++ -g -o simpletest simpletest.cpp electroweakino.o \
`root-config --cflags --glibs` -lblas -llapack -lgfortran

rm lprint.inc

echo 'findmasses executable maker seems to have completed OK'

exit
