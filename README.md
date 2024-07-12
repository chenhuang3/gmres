This is a conversion of the MATLAB GMRES from netlib https://www.netlib.org/templates/matlab/gmres.m

Some main modifications are:

1. A*x is replaced by an external user-defined subroutine (matvec).
   
3. Norm and dot product of vectors can take a metric, that is, sum(v1 * v2) is replaced by sum(v1 * v2) * metric. If no metric is needed, simply set metric=1.d0.

An example for using the code is given in test_gmres.f90
