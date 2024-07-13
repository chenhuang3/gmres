This is a conversion of the MATLAB GMRES for solving A*x=b from netlib https://www.netlib.org/templates/matlab/gmres.m

A main modification is that A * x is replaced by an external user-defined subroutine (matvec).

An example of using the code is given in test_gmres.f90. The test can be complied using ./compile.sh. Lapack is needed.

For any comment/suggestions, please email to chuang3@fsu.edu  
