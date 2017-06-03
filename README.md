# kdrunoff

Use the kdtree algorithm to remap river runoff from land to ocean. Fortran and Python. 

## To run the tests

```
gfortran -fdefault-real-8 kdtree2/src-f90/kdtree2.f90 kdrunoff.F90 test.F90 -o test.exe -lnetcdf -lnetcdff
./test.exe
```
