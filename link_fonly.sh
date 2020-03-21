#!/bin/sh
# 
# This is the Fortran only version for a 
# specific (M1, M2, mu, tanb) point.
#
# Compile the main program, ewdriver.f90 
# with the computational subroutine electroweakino.f90 
# and link with LAPACK and BLAS libraries
#
# This could potentially use command line arguments or 
# an input file to specify the points, or an include file.
# Latter is the most straightforward for now.
#
#                Graham W. Wilson      Reviewed again 06-Sep-2019
#
ln -s lprint_true.inc lprint.inc

ln -s ewdriver_Method2_Point0.dat ewdriver_inputs.dat
#ln -s glpA.dat ewdriver_inputs.dat
#ln -s glpD.dat ewdriver_inputs.dat
#ln -s test.dat ewdriver_inputs.dat

gfortran -c electroweakino.f90
gfortran -o ewdriver ewdriver.f90 electroweakino.o -lblas -llapack
rm electroweakino.o   #not needed now
rm lprint.inc

rm ewdriver_inputs.dat

exit
